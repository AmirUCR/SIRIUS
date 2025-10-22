#pragma once
#include "tables.hpp"

#include <set>
#include <string>
#include <vector>
#include <unordered_map>

namespace sirius
{
    class SIRIUSInstance
    {
    public:
        int    n = 0;
        int    decidable_protein_length = 0;
        int    warm_start_largest_fragment_length = 0;
        
        double hard_rscu_threshold = 0.0;
        double soft_rscu_threshold = 0.0;
        double gc_ending_rscu_threshold = 0.0;
        double rscu_alpha = 0.0;
        double max_low_rscu_ratio = 0.0;

        bool   hard_filter_by_rscu = false;
        bool   soft_filter_by_rscu = false;
        
        SIRIUSTables tables;
        
        std::string dna_with_holes;
        std::string decidable_protein;
        std::string init_target_protein;
        std::string codon_usage_path;
        std::string output_folder;
        
        std::vector<std::vector<std::string>> warm_start_solution_seqs;

        SIRIUSInstance() = default;

        SIRIUSInstance(
            int n,
            std::string init_target_protein,
            const SIRIUSTables& tables,
            std::string codon_usage_path,
            std::string output_folder,
            double hard_rscu_threshold,
            double soft_rscu_threshold,
            double gc_ending_rscu_threshold,
            double rscu_alpha,
            double max_low_rscu_ratio
        );

        void recompute_flags();
        
        void input_to_uppercase();
        
        std::vector<std::vector<std::string>> process_warm_start_output(
            const std::string &python_output,
            const std::unordered_map<char, std::string> &codon_table,
            const std::unordered_map<std::string, char> &translate_codon_table,
            const std::set<char> &skip_aa_table,
            int codon_length
        );

        int find_max_fragment_length_from_seqs(
            const std::vector<std::vector<std::string>> &seqs);
    };
}