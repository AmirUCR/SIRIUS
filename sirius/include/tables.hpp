#pragma once
#include <set>
#include <string>
#include <vector>
#include <unordered_map>

namespace sirius
{
    class SIRIUSTables
    {
    public:
        SIRIUSTables();

        std::set<char> skip_aa;
        std::unordered_map<char, int> reduced_codon_lengths_table;
        std::unordered_map<char, std::string> invariant_codon_table;
        std::unordered_map<std::string, char> translate_codon_table;
        std::unordered_map<char, std::vector<std::string>> reduced_codon_table;
        std::unordered_map<char, std::unordered_map<std::string, double>> rscu_map;

        void build_rscu_map_from_csv(const std::string &filename);
    };
}