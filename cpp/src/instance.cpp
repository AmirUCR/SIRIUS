#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "utils.hpp"
#include "tables.hpp"
#include "logging.hpp"
#include "instance.hpp"

namespace sirius 
{
    SIRIUSInstance::SIRIUSInstance(
        int n,
        std::string init_target_protein,
        const SIRIUSTables& tables,
        std::string codon_usage_path,
        std::string output_folder,
        double hard_thr,
        double soft_thr,
        double gc_end_thr,
        double alpha,
        double max_ratio)
    : n(n)
    , decidable_protein_length(0)
    , warm_start_largest_fragment_length(0)
    , hard_rscu_threshold(hard_thr)
    , soft_rscu_threshold(soft_thr)
    , gc_ending_rscu_threshold(gc_end_thr)
    , rscu_alpha(alpha)
    , max_low_rscu_ratio(max_ratio)
    , hard_filter_by_rscu(false)
    , soft_filter_by_rscu(false)
    , tables(std::move(tables))
    , dna_with_holes("")
    , decidable_protein("")
    , init_target_protein(std::move(init_target_protein))
    , codon_usage_path(std::move(codon_usage_path))
    , output_folder(std::move(output_folder))
    , warm_start_solution_seqs()
    {
        recompute_flags();
        input_to_uppercase();

        double threshold_genediversifier = this->hard_rscu_threshold;
        if (this->soft_filter_by_rscu)
        {
            threshold_genediversifier = this->soft_rscu_threshold;
        }

        std::string py_output = utils::call_gene_diversifier(
            this->init_target_protein, 
            this->n,
            this->codon_usage_path,
            threshold_genediversifier,
            gc_ending_rscu_threshold);

        this->warm_start_solution_seqs = process_warm_start_output(
            py_output,
            this->tables.invariant_codon_table,
            this->tables.translate_codon_table,
            this->tables.skip_aa,
            3);

        this->warm_start_largest_fragment_length = find_max_fragment_length_from_seqs(
            this->warm_start_solution_seqs);

        logg.print_info_newline("Warm-start from fragment of length " +
            std::to_string(this->warm_start_largest_fragment_length));
    }

    void SIRIUSInstance::input_to_uppercase() {
        std::transform(
            this->init_target_protein.begin(),
            this->init_target_protein.end(),
            this->init_target_protein.begin(), ::toupper);
            
        for (const char& c : this->init_target_protein) {
            this->dna_with_holes += this->tables.invariant_codon_table[c];

            if (this->tables.skip_aa.find(c) == this->tables.skip_aa.end())
            {   // if c not in skip_aa
                this->decidable_protein += c;
                this->decidable_protein_length += 1;
            }
        }
    }

    void SIRIUSInstance::recompute_flags() {
        hard_filter_by_rscu = (hard_rscu_threshold > 0.0);
        soft_filter_by_rscu = (soft_rscu_threshold > 0.0);
    }

    std::vector<std::vector<std::string>>
    SIRIUSInstance::process_warm_start_output(
        const std::string &python_output,
        const std::unordered_map<char, std::string> &codon_table,
        const std::unordered_map<std::string, char> &translate_codon_table,
        const std::set<char> &skip_aa_table,
        int codon_length)
    {
        std::istringstream in(python_output);
        std::string line;
        std::vector<std::vector<std::string>> seqs;

        while (std::getline(in, line))
        {
            // Trim possible trailing carriage returns or spaces
            if (!line.empty() && (line.back() == '\r' || line.back() == '\n')) {
                line.pop_back();
            }
            if (line.empty()) {
                continue;
            }

            std::vector<std::string> variable_bases_of_codons;

            for (size_t i = 0; i + codon_length - 1 < line.size(); i += codon_length)
            {
                std::string codon = line.substr(i, codon_length);
                char aa = translate_codon_table.at(codon);

                if (skip_aa_table.find(aa) != skip_aa_table.end())
                {
                    continue;
                }

                std::string bases = codon_table.at(aa);
                std::string variable_bases;
                for (size_t base_i = 0; base_i < bases.size(); ++base_i)
                {
                    if (bases[base_i] == '_')
                    {
                        variable_bases += codon[base_i];
                    }
                }

                if (!variable_bases.empty())
                {
                    variable_bases_of_codons.push_back(variable_bases);
                }
            }

            seqs.push_back(variable_bases_of_codons);
        }

        return seqs;
    }

    int SIRIUSInstance::find_max_fragment_length_from_seqs(
        const std::vector<std::vector<std::string>> &seqs)
    {
        int max_length = 0;

        // Flatten each warm-start sequence (vector<string> --> string)
        std::vector<std::string> flat_seqs;
        for (const auto &seq : seqs)
        {
            std::string flat;
            for (const std::string &token : seq)
            {
                flat += token;
            }
            flat_seqs.push_back(flat);
        }

        // Compare every unique pair of sequences
        for (size_t i = 0; i < flat_seqs.size(); ++i)
        {
            for (size_t j = i + 1; j < flat_seqs.size(); ++j)
            {
                const std::string &s1 = flat_seqs[i];
                const std::string &s2 = flat_seqs[j];

                int current = 0;
                int local_max = 0;
                for (size_t k = 0; k < std::min(s1.size(), s2.size()); ++k)
                {
                    if (s1[k] == s2[k])
                    {
                        ++current;
                        local_max = std::max(local_max, current);
                    }
                    else
                    {
                        current = 0;
                    }
                }

                max_length = std::max(max_length, local_max);
            }
        }

        return max_length;
    }
}