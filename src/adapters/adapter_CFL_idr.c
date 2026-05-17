#include "adapter_CFL_idr.h"

#include "parser.h"

#define TRY(GrB_method)                                                                                                \
    {                                                                                                                  \
        GrB_Info LG_GrB_Info = GrB_method;                                                                             \
        if (LG_GrB_Info < GrB_SUCCESS) {                                                                               \
            fprintf(stderr, "LAGraph failure (file %s, line %d): (%d) \n", __FILE__, __LINE__, LG_GrB_Info);           \
            return (LG_GrB_Info);                                                                                      \
        }                                                                                                              \
    }

// inner state of the adapter
typedef struct {
    IdrGraph working_graph;
    IdrGraph parsed_graph;
    bool working_graph_owned;
    GrB_Matrix under_result;
    GrB_Matrix over_result;
    GrB_Matrix on_demand_result;
    IdrPrepareData opts;
} state_t;

static state_t state;

static GrB_Info build_idr_graph_from_parser(IdrGraph *out, ParserResult pr) {
    GraphMatrices gm = get_grb_matrices_from_graph(pr.graph, &pr.symbols);

    GrB_Info info = get_idr_graph(out, &gm, &pr.symbols, pr.node_count, &DefaultTerminalFormat);

    if (gm.matrices != NULL) {
        for (size_t i = 0; i < gm.count; i++) {
            if (gm.matrices[i] != NULL) {
                GrB_Matrix_free(&gm.matrices[i]);
            }
        }
        free((void *)gm.matrices);
    }
    free(gm.matrix_symbols);

    return info;
}

// initialize LAGraph/GraphBLAS
static GrB_Info adapter_CFL_idr_setup() {
    char msg[LAGRAPH_MSG_LEN];
    TRY(LAGraph_Init(msg));
    return GrB_SUCCESS;
}

// inner structure to pass data to prepare function
typedef IdrPrepareData PrepareData;

// prepare the adapter for use with the given parser result
// this may include converting the graph into a format suitable for the algorithm
//
// this modifies the inner state of the adapter
//
// adapter_CFL_idr_prepare should be called just once for each config
static GrB_Info adapter_CFL_idr_prepare(ParserResult parser_result, void *prepare_data) {
    PrepareData *data = (PrepareData *)prepare_data;
    if (data) {
        state.opts = *data;
    }

    TRY(build_idr_graph_from_parser(&state.parsed_graph, parser_result));

    if (state.opts.valueflow) {
        TRY(idr_remove_valueflow_unreachable(&state.working_graph, &state.parsed_graph));
        idr_graph_free(&state.parsed_graph);
        state.working_graph_owned = true;
    } else {
        state.working_graph = state.parsed_graph;
        state.working_graph_owned = false;
    }

    return GrB_SUCCESS;
}

// initialize output matrices
//
// this should be called before each run of the algorithm
static GrB_Info adapter_CFL_idr_init_outputs() {
    state.under_result = NULL;
    state.over_result = NULL;
    state.on_demand_result = NULL;
    return GrB_SUCCESS;
}

// run the algorithm
//
// this should be called after adapter_CFL_idr_init_outputs
static GrB_Info adapter_CFL_idr_run() {
    TRY(idr_get_under_approx(&state.under_result, &state.working_graph, state.opts.valueflow));
    GrB_Matrix_wait(state.under_result, GrB_MATERIALIZE);

    TRY(idr_get_over_approx(&state.over_result, &state.working_graph, state.opts.grammar_type, state.under_result,
                            state.opts.valueflow, /*filter_empty=*/true));
    GrB_Matrix_wait(state.over_result, GrB_MATERIALIZE);

    if (state.opts.on_demand || state.opts.parity_d) {
        TRY(idr_get_on_demand(&state.on_demand_result, &state.working_graph, state.under_result, state.over_result,
                              state.opts.parity_d, state.opts.valueflow, /*filter_empty=*/true));
        GrB_Matrix_wait(state.on_demand_result, GrB_MATERIALIZE);
    }

    return GrB_SUCCESS;
}

// check if the result of the algorithm is valid
static ResultType adapter_CFL_idr_is_result_valid(size_t valid_result) {
    GrB_Matrix target = state.on_demand_result ? state.on_demand_result : state.over_result;
    GrB_Index nvals = 0;
    if (target) {
        GrB_Matrix_nvals(&nvals, target);
    }
    return (size_t)nvals == valid_result ? RESULT_OK : RESULT_ERROR;
}

// get the result of the algorithm
static size_t adapter_CFL_idr_get_result() {
    GrB_Matrix target = state.on_demand_result ? state.on_demand_result : state.over_result;
    GrB_Index nvals = 0;
    if (target) {
        GrB_Matrix_nvals(&nvals, target);
    }
    return (size_t)nvals;
}

// free output matrices
//
// this should be called after each run of the algorithm
static GrB_Info adapter_CFL_idr_free_outputs() {
    GrB_Matrix_free(&state.under_result);
    GrB_Matrix_free(&state.over_result);
    GrB_Matrix_free(&state.on_demand_result);
    return GrB_SUCCESS;
}

// free all resources allocated by the adapter except outputs
//
// this should be called after all runs of the algorithm for the given config
static GrB_Info adapter_CFL_idr_cleanup() {
    if (state.working_graph_owned) {
        idr_graph_free(&state.working_graph);
    } else {
        idr_graph_free(&state.parsed_graph);
    }
    state.working_graph_owned = false;
    return GrB_SUCCESS;
}

// free LAGraph/GraphBLAS resources
static GrB_Info adapter_CFL_idr_teardown() {
    char msg[LAGRAPH_MSG_LEN];
    TRY(LAGraph_Finalize(msg));
    return GrB_SUCCESS;
}

// get the methods of the adapter
AdapterMethods adapter_CFL_idr_get_methods() {
    return (AdapterMethods){
        .setup = adapter_CFL_idr_setup,
        .teardown = adapter_CFL_idr_teardown,
        .init_outputs = adapter_CFL_idr_init_outputs,
        .free_outputs = adapter_CFL_idr_free_outputs,
        .run = adapter_CFL_idr_run,
        .prepare = adapter_CFL_idr_prepare,
        .cleanup = adapter_CFL_idr_cleanup,
        .is_result_valid = adapter_CFL_idr_is_result_valid,
        .get_result = adapter_CFL_idr_get_result,
    };
}
