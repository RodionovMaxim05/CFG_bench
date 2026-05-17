#pragma once

#include "adapter.h"
#include "cfl_idr_bench.h"

typedef struct IdrPrepareData {
    bool valueflow;
    bool on_demand;
    bool parity_d;
    IdrGrammarType grammar_type;
} IdrPrepareData;

AdapterMethods adapter_CFL_idr_get_methods(void);
