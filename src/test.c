#include "adapter_CFL.h"
#include "adapter_CFL_CFPQ_RSM.h"
#include "adapter_CFL_adv.h"
#include "adapter_CFL_all_path.h"
#include "adapter_CFL_all_path_adv.h"
#include "adapter_CFL_idr.h"
#include "adapter_CFL_multsrc.h"
#include "adapter_CFL_single_path.h"
#include "memory.h"
#include "parser.h"
#include "result_manager.h"
#include "stats.h"
#include <GraphBLAS.h>
#include <LAGraph.h>
#include <LAGraphX.h>
#include <getopt.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define run_algorithm()                                                                                                \
    LAGraph_CFL_reachability_adv(outputs, adj_matrices, symbols_amount, grammar.rules, grammar.rules_count, msg,       \
                                 optimizations)

#define check_error(error)                                                                                             \
    {                                                                                                                  \
        retval = run_algorithm();                                                                                      \
        TEST_CHECK(retval == error);                                                                                   \
        TEST_MSG("retval = %d (%s)", retval, msg);                                                                     \
    }

#define check_result(result)                                                                                           \
    {                                                                                                                  \
        char *expected = output_to_str(0);                                                                             \
        TEST_CHECK(strcmp(result, expected) == 0);                                                                     \
        TEST_MSG("Wrong result. Actual: %s", expected);                                                                \
    }

#define TRY(GrB_method)                                                                                                \
    {                                                                                                                  \
        GrB_Info LG_GrB_Info = GrB_method;                                                                             \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                                               \
            fprintf(stderr, "LAGraph failure (file %s, line %d): (%d) \n", __FILE__, __LINE__, LG_GrB_Info);           \
            return (LG_GrB_Info);                                                                                      \
        }                                                                                                              \
    }

GrB_Matrix *adj_matrices = NULL;
GrB_Matrix *outputs = NULL;
grammar_t grammar = {0, 0, NULL};
char msg[LAGRAPH_MSG_LEN];
size_t symbols_amount = 0;

void print_rules(Grammar grammar, SymbolList list) {
    for (size_t i = 0; i < grammar.rules_count; i++) {
        Rule rule = grammar.rules[i];
        if (rule.first != -1) {
            printf("%s ->", list.symbols[rule.first].label);
        }
        if (rule.second != -1) {
            printf(" %s", list.symbols[rule.second].label);
        }
        if (rule.third != -1) {
            printf(" %s", list.symbols[rule.third].label);
        }
        printf("\n");
    }
}

void print_list(SymbolList list, size_t *map) {
    if (map == NULL) {
        for (size_t i = 0; i < list.count; i++) {
            Symbol sym = list.symbols[i];
            printf("[%2ld] %s [%s] %s\n", i, sym.label, sym.is_nonterm ? "N" : "T", sym.is_indexed ? "[I]" : "");
        }
    } else {
        for (size_t i = 0; i < list.count; i++) {
            Symbol sym = list.symbols[i];
            printf("[%2ld] %s [%s] %s\n", map[i], sym.label, sym.is_nonterm ? "N" : "T", sym.is_indexed ? "[I]" : "");
        }
    }
}

#define RESET "\033[0m"
#define BLACK "\033[30m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define WHITE "\033[37m"

// Use your custom configuration for the benchmark (default is the xz.g graph
// and vf.cnf grammar)
#define configs_macro configs_java

#define OPT_EMPTY (1 << 0)
#define OPT_FORMAT (1 << 1)
#define OPT_LAZY (1 << 2)
#define OPT_BLOCK (1 << 3)

enum { HOT_OPTION = 1000, FULL_TIME_OPTION };

typedef struct {
    AdapterMethods adapter;
    char *algo;
    char *input_config;
    size_t rounds_count;
    bool is_test;
    bool is_hot_enabled;
    int8_t optimizations;
    IdrPrepareData idr_opts;
} BenchArgs;

static void print_usage(const char *program_name) {
    fprintf(stderr,
            "Usage: %s -c <config file> [options]\n"
            "\n"
            "Required:\n"
            "  -c <config file>  Path to benchmark config file\n"
            "\n"
            "Benchmark options:\n"
            "  -r <rounds>       Number of benchmark rounds (default: 10)\n"
            "  --hot             Enable HOT launch (warm-up run before measurements)\n"
            "  --full-time       Measure full time including parsing and prepare\n"
            "  -a <algorithm>    Algorithm to use "
            "(default: CFL_adv; options: CFL_adv, CFL, CFL_single_path, CFL_all_path, CFL_all_path_adv, CFL_CFPQ_RSM, "
            "CFL_multsrc, CFL_idr)\n"
            "\n"
            "Optimization flags:\n"
            "  -e                Enable empty optimization\n"
            "  -f                Enable format optimization\n"
            "  -l                Enable lazy optimization\n"
            "  -b                Enable block optimization\n"
            "\n"
            "Other:\n"
            "  -t                Enable test mode\n"
            "  -h                Print this help message\n"
            "\n"
            "Example:\n"
            "  %s -c configs/configs_my.csv -r 10 --hot\n",
            program_name, program_name);
}

static void resolve_prepare_data(const BenchArgs *args, void **prepare_data, CFL_adv_PrepareData *cfl_adv_data) {
    if (strcmp(args->algo, "CFL_idr") == 0) {
        *prepare_data = (void *)&args->idr_opts;
    } else {
        cfl_adv_data->optimizations = args->optimizations;
        *prepare_data = cfl_adv_data;
    }
}

static void resolve_idr_opts(IdrPrepareData *idr_opts, const char *algo, config_row config) {
    if (strcmp(algo, "CFL_idr") != 0 || config.idr_grammar == NULL) {
        return;
    }

    idr_opts->parity_d = false;
    idr_opts->on_demand = false;

    if (strcmp(config.idr_grammar, "parity") == 0)
        idr_opts->grammar_type = IDR_PARITY;
    else if (strcmp(config.idr_grammar, "parity2") == 0)
        idr_opts->grammar_type = IDR_PARITY2;
    else if (strcmp(config.idr_grammar, "se") == 0)
        idr_opts->grammar_type = IDR_SE;
    else if (strcmp(config.idr_grammar, "project") == 0)
        idr_opts->grammar_type = IDR_PROJECT;
    else if (strcmp(config.idr_grammar, "exclude") == 0)
        idr_opts->grammar_type = IDR_EXCLUDE;
    else if (strcmp(config.idr_grammar, "all") == 0)
        idr_opts->grammar_type = IDR_ALL;
    else if (strcmp(config.idr_grammar, "parityD") == 0) {
        idr_opts->grammar_type = IDR_PARITY;
        idr_opts->parity_d = true;
    } else if (strcmp(config.idr_grammar, "on-demand") == 0) {
        idr_opts->grammar_type = IDR_ALL;
        idr_opts->on_demand = true;
    } else {
        fprintf(stderr, "Unknown idr grammar: %s\n", config.idr_grammar);
        exit(EXIT_FAILURE);
    }
}

static void print_test_result(size_t result, GrB_Info retval, config_row config, AdapterMethods adapter,
                              double elapsed) {
    ResultType result_type = adapter.is_result_valid(config.valid_result);
    char status[256];
    switch (result_type) {
    case RESULT_OK:
        snprintf(status, sizeof(status), GREEN "[OK]" RESET);
        break;
    case RESULT_ERROR:
        snprintf(status, sizeof(status), RED "[Wrong] (Result must be %ld)" RESET, config.valid_result);
        break;
    case RESULT_UNKNOWN:
        snprintf(status, sizeof(status), YELLOW "[Unknown]" RESET);
        break;
    default:
        fprintf(stderr, "Unknown result type: %d\n", result_type);
        abort();
    }
    printf("\tResult: %ld (Return code: %d) %s (%.4f sec)", result, retval, status, elapsed);
}

typedef enum { BENCH_ALGO_ONLY, BENCH_FULL_TIME } BenchMode;

static GrB_Info bench(BenchArgs *args, config_row *configs, size_t configs_count, BenchMode mode) {
    GrB_Info retval = GrB_SUCCESS;

    double setup_time = LAGraph_WallClockTime();
    TRY(args->adapter.setup());
    setup_time = LAGraph_WallClockTime() - setup_time;

    for (size_t i = 0; i < configs_count; i++) {
        double *elapsed = calloc(args->rounds_count, sizeof(double));
        ssize_t *memory_kb = calloc(args->rounds_count, sizeof(ssize_t));
        if (!elapsed || !memory_kb) {
            free(elapsed);
            free(memory_kb);
            return GrB_OUT_OF_MEMORY;
        }

        config_row config = configs[i];
        printf("CONFIG: grammar: %s, graph: %s\n", config.idr_grammar ? config.idr_grammar : config.grammar,
               config.graph);
        fflush(stdout);

        resolve_idr_opts(&args->idr_opts, args->algo, config);

        void *prepare_data = NULL;
        CFL_adv_PrepareData cfl_adv_data = {0};
        resolve_prepare_data(args, &prepare_data, &cfl_adv_data);

        // In algo-only mode, parse once before the loop
        ParserResult parser_result_outer = {0};
        if (mode == BENCH_ALGO_ONLY) {
            parser_result_outer = parser(config);
            args->adapter.prepare(parser_result_outer, prepare_data);
        }

        bool is_hot = args->is_hot_enabled;

        size_t result = 0;
        for (size_t j = 0; j < args->rounds_count; j++) {
            if (mode == BENCH_ALGO_ONLY) {
                TRY(args->adapter.init_outputs());
            }

            // in some cases free don't change memory usage, so we need to reset it manually
            malloc_trim(0);
            if (mem_peak_reset() != 0) {
                fprintf(stderr, "Failed to reset memory peak\n");
                exit(EXIT_FAILURE);
            }

            double t_start = LAGraph_WallClockTime();

            if (mode == BENCH_FULL_TIME) {
                ParserResult pr = parser(config);
                args->adapter.prepare(pr, prepare_data);
                TRY(args->adapter.init_outputs());
            }

#ifndef CI
            retval = args->adapter.run();
#endif
            double t_end = LAGraph_WallClockTime();
            memory_kb[j] = mem_get_peak_kb();

            double round_time = (t_end - t_start);
            if (mode == BENCH_FULL_TIME) {
                round_time += setup_time;
            }
            elapsed[j] = round_time;

            if (args->is_test) {
                result = args->adapter.get_result();
                print_test_result(result, retval, config, args->adapter, round_time);
                if (retval != 0) {
                    printf("\t(MSG: %s)", msg);
                }
                printf("\n");

                TRY(args->adapter.free_outputs());
                if (mode == BENCH_FULL_TIME) {
                    args->adapter.cleanup();
                }
                break;
            }

            if (is_hot) {
                is_hot = false;
                j--;
                TRY(args->adapter.free_outputs());
                if (mode == BENCH_FULL_TIME) {
                    args->adapter.cleanup();
                }
                continue;
            }

            printf("\t%.3f s %.3f MB\n", round_time, (double)memory_kb[j] / 1024);
            fflush(stdout);

            result = args->adapter.get_result();
            TRY(args->adapter.free_outputs());

            if (mode == BENCH_FULL_TIME) {
                args->adapter.cleanup();
            }

            save_result(args->algo, config.grammar, config.graph, result, memory_kb[j], (size_t)(round_time * 1000));
            // in some cases free don't change memory usage, so we need to reset it manually
            malloc_trim(0);
        }
        printf("\n");

        if (mode == BENCH_ALGO_ONLY && !args->is_test) {
            TRY(args->adapter.cleanup());
        }

        if (!args->is_test) {
            StatResult time_stats = calc_stats(elapsed, args->rounds_count);

            double *memory_mb = malloc(args->rounds_count * sizeof(double));
            if (!memory_mb) {
                fprintf(stderr, "Failed to allocate memory for stats\n");
                free(elapsed);
                free(memory_kb);
                return GrB_OUT_OF_MEMORY;
            }
            for (size_t j = 0; j < args->rounds_count; j++) {
                memory_mb[j] = (double)(memory_kb[j]) / 1024;
            }
            StatResult mem_stats = calc_stats(memory_mb, args->rounds_count);
            free(memory_mb);

            printf("\tTime: %.6f ± %.6f s "
                   "Memory: %.3f ± %.3f MB "
                   "Result: %ld (return code %d) (%s)\n\n",
                   time_stats.mean, time_stats.error, mem_stats.mean, mem_stats.error, result, retval, msg);
        }

        free(elapsed);
        free(memory_kb);
        fflush(stdout);
    }

    return GrB_SUCCESS;
}

int main(int argc, char **argv) {
    int opt;
    bool is_config = false;
    bool is_algo_chosen = false;
    bool is_full_time = false;
    char *input_config = NULL;

    BenchArgs args = {
        .rounds_count = 10,
        .is_test = false,
        .is_hot_enabled = false,
        .optimizations = 0,
        .idr_opts =
            {
                .valueflow = false,
                .on_demand = false,
                .parity_d = false,
                .grammar_type = IDR_DEFAULT,
            },
    };

    static struct option long_options[] = {
        {"hot", no_argument, 0, HOT_OPTION}, {"full-time", no_argument, 0, FULL_TIME_OPTION}, {0, 0, 0, 0}};

    while ((opt = getopt_long(argc, argv, "eflbthr:c:a:", long_options, NULL)) != -1) {
        switch (opt) {
        case 'e':
            args.optimizations |= OPT_EMPTY;
            break;
        case 'f':
            args.optimizations |= OPT_FORMAT;
            break;
        case 'l':
            args.optimizations |= OPT_LAZY;
            break;
        case 'b':
            args.optimizations |= OPT_BLOCK;
            break;
        case 'h':
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
        case HOT_OPTION:
            args.is_hot_enabled = true;
            break;
        case FULL_TIME_OPTION:
            is_full_time = true;
            break;
        case 't':
            args.is_test = true;
            break;
        case 'r':
            args.rounds_count = strtoul(optarg, NULL, 10);
            if (args.rounds_count == 0) {
                fprintf(stderr, "Rounds count must be greater than 0\n");
                exit(EXIT_FAILURE);
            }
            printf("Choosen rounds count: %zu\n", args.rounds_count);
            break;
        case 'c':
            is_config = true;
            input_config = optarg;
            printf("Choosen config: %s\n", input_config);
            break;
        case 'a':
            is_algo_chosen = true;
            args.algo = optarg;
            printf("Choosen algorithm: %s\n", args.algo);

            if (strcmp(args.algo, "CFL_adv") == 0)
                args.adapter = adapter_CFL_adv_get_methods();
            else if (strcmp(args.algo, "CFL") == 0)
                args.adapter = adapter_CFL_get_methods();
            else if (strcmp(args.algo, "CFL_single_path") == 0)
                args.adapter = adapter_CFL_single_path_get_methods();
            else if (strcmp(args.algo, "CFL_all_path") == 0)
                args.adapter = adapter_CFL_all_paths_get_methods();
            else if (strcmp(args.algo, "CFL_CFPQ_RSM") == 0)
                args.adapter = adapter_CFL_CFPQ_RSM_get_methods();
            else if (strcmp(args.algo, "CFL_multsrc") == 0)
                args.adapter = adapter_CFL_multsrc_get_methods();
            else if (strcmp(args.algo, "CFL_all_path_adv") == 0)
                args.adapter = adapter_CFL_all_path_adv_get_methods();
            else if (strcmp(args.algo, "CFL_idr") == 0)
                args.adapter = adapter_CFL_idr_get_methods();
            else {
                fprintf(stderr, "Unknown algorithm: %s\n", args.algo);
                exit(EXIT_FAILURE);
            }
            break;
        default:
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (!is_algo_chosen) {
        args.adapter = adapter_CFL_adv_get_methods();
        args.algo = "CFL_adv";
        printf("No algorithm chosen, using CFL_adv by default\n");
    }

    if (!is_config) {
        fprintf(stderr, "Need to choose config by flag -c [config file]\n");
        exit(EXIT_FAILURE);
    }

    size_t configs_count = 0;
    config_row *configs = calloc(1000, sizeof(config_row));
    char *config_text;
    get_configs_from_file(input_config, &configs_count, configs, &config_text);

    printf("Start bench (%s)\n", is_full_time ? "full-time mode" : "algorithm-only mode");
    fflush(stdout);

    GrB_Info retval = bench(&args, configs, configs_count, is_full_time ? BENCH_FULL_TIME : BENCH_ALGO_ONLY);

    free(configs);
    free(config_text);
    TRY(args.adapter.teardown());
    return retval == GrB_SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
