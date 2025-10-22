#pragma once
#include <cstdint>
#include "ortools/sat/sat_parameters.pb.h"

namespace sirius
{
    class SIRIUSConfig 
    {
    public:
        operations_research::sat::SatParameters parameters;

        SIRIUSConfig();  // default

        SIRIUSConfig(
            bool log_search_progress,
            unsigned int num_workers,
            unsigned int linearization_level,
            double relative_gap_limit,
            int max_time_in_seconds,
            bool fix_variables_to_their_hinted_value = false);
    };
} // namespace sirius
