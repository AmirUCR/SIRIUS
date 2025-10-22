// https://github.com/google/or-tools/blob/stable/ortools/sat/sat_parameters.proto
#include "config.hpp"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/sat_parameters.pb.h"

namespace sirius 
{
    SIRIUSConfig::SIRIUSConfig() = default;

    SIRIUSConfig::SIRIUSConfig(
        bool log_search_progress,
        unsigned int num_workers,
        unsigned int linearization_level,
        double relative_gap_limit,
        int max_time_in_seconds,
        bool fix_variables_to_their_hinted_value)
    {
        parameters.set_num_workers(static_cast<int>(num_workers));
        parameters.set_log_search_progress(log_search_progress);
        parameters.set_linearization_level(static_cast<int>(linearization_level));
        parameters.set_relative_gap_limit(relative_gap_limit);
        parameters.set_max_time_in_seconds(max_time_in_seconds);
        parameters.set_fix_variables_to_their_hinted_value(fix_variables_to_their_hinted_value);
    }
} // namespace sirius
