#pragma once
#include <utility>

#include "config.hpp"
#include "tables.hpp"
#include "instance.hpp"

namespace sirius 
{
    class SIRIUSInterface 
    {
    public:
        std::pair<SIRIUSInstance, SIRIUSConfig>
            gather_inputs_from_flags(int argc, char* argv[], const SIRIUSTables& tables);

        std::pair<SIRIUSInstance, SIRIUSConfig>
            gather_inputs_interactively(const SIRIUSTables& tables);
    };
}