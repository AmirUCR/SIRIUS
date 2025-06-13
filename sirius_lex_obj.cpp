
#include <map>
#include <queue>
#include <deque>
#include <tuple>
#include <memory>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include <unordered_map>
#include <unordered_set>

#include "ortools/util/sigint.h"
#include "ortools/base/logging.h"
#include "ortools/base/init_google.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/sat_parameters.pb.h"

#define RESET "\033[0m"
#define RED "\033[31m"
#define BLUE "\033[34;1m"
#define ORANGE "\033[38;5;208m"

std::string elapsed_since_start();

using VectorBoolVariables = std::vector<operations_research::sat::BoolVar>;
using VectorLinearExprs = std::vector<operations_research::sat::LinearExpr>;
using VectorColumnYVars = std::vector<VectorLinearExprs>;
using VectorColumnZVars = std::vector<VectorBoolVariables>;
using VectorCodonBoolVariables = std::vector<VectorBoolVariables>;
using VectorAminoacidCodonBoolVariables = std::vector<VectorCodonBoolVariables>;
using VectorSequenceAminoacidCodonBoolVariables = std::vector<VectorAminoacidCodonBoolVariables>;

// Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;

struct SIRIUSTables
{
    // std::unordered_map<char, std::string> invariant_codon_table = {
    //     {'A', "GC_"}, // {"GCT", "GCC", "GCA", "GCG"}
    //     {'C', "TG_"}, // {"TGT", "TGC"}
    //     {'D', "GA_"}, // {"GAT", "GAC"}
    //     {'E', "GA_"}, // {"GAA", "GAG"}
    //     {'F', "TT_"}, // {"TTT", "TTC"}
    //     {'G', "GG_"}, // {"GGT", "GGC", "GGA", "GGG"}
    //     {'H', "CA_"}, // {"CAT", "CAC"}
    //     {'I', "AT_"}, // {"ATT", "ATC", "ATA"}
    //     {'K', "AA_"}, // {"AAA", "AAG"}
    //     {'L', "_T_"}, // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    //     {'M', "ATG"},
    //     {'N', "AA_"}, // {"AAT", "AAC"}
    //     {'P', "CC_"}, // {"CCT", "CCC", "CCA", "CCG"}
    //     {'Q', "CA_"}, // {"CAA", "CAG"}
    //     {'R', "_G_"}, // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
    //     {'S', "___"}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
    //     {'T', "AC_"}, // {"ACT", "ACC", "ACA", "ACG"}
    //     {'V', "GT_"}, // {"GTT", "GTC", "GTA", "GTG"}
    //     {'W', "TGG"}, // TGG
    //     {'Y', "TA_"}, // {"TAT", "TAC"}
    //     {'*', "T__"}  // {"TAA", "TAG", "TGA"}
    // };

    std::unordered_map<char, std::string> invariant_codon_table = {
        {'A', "___"}, // {"GCT", "GCC", "GCA", "GCG"}
        {'C', "___"}, // {"TGT", "TGC"}
        {'D', "___"}, // {"GAT", "GAC"}
        {'E', "___"}, // {"GAA", "GAG"}
        {'F', "___"}, // {"TTT", "TTC"}
        {'G', "___"}, // {"GGT", "GGC", "GGA", "GGG"}
        {'H', "___"}, // {"CAT", "CAC"}
        {'I', "___"}, // {"ATT", "ATC", "ATA"}
        {'K', "___"}, // {"AAA", "AAG"}
        {'L', "___"}, // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', "___"},
        {'N', "___"}, // {"AAT", "AAC"}
        {'P', "___"}, // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', "___"}, // {"CAA", "CAG"}
        {'R', "___"}, // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', "___"}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', "___"}, // {"ACT", "ACC", "ACA", "ACG"}
        {'V', "___"}, // {"GTT", "GTC", "GTA", "GTG"}
        {'W', "___"}, // TGG
        {'Y', "___"}, // {"TAT", "TAC"}
        {'*', "___"}  // {"TAA", "TAG", "TGA"}
    };

    // std::set<char> skip_aa = {'M', 'W'};
    std::set<char> skip_aa = {};

    // std::unordered_map<char, std::vector<std::string>> reduced_codon_table = {
    //     {'A', {"A", "C", "G", "T"}},                       // {"GCT", "GCC", "GCA", "GCG"}
    //     {'C', {"C", "T"}},                                 // {"TGT", "TGC"}
    //     {'D', {"C", "T"}},                                 // {"GAT", "GAC"}
    //     {'E', {"A", "G"}},                                 // {"GAA", "GAG"}
    //     {'F', {"C", "T"}},                                 // {"TTT", "TTC"}
    //     {'G', {"A", "C", "G", "T"}},                     // {"GGT", "GGC", "GGA", "GGG"}
    //     {'H', {"C", "T"}},                                 // {"CAT", "CAC"}
    //     {'I', {"A", "C", "T"}},                            // {"ATT", "ATC", "ATA"}
    //     {'K', {"A", "G"}},                                 // {"AAA", "AAG"}
    //     {'L', {"CA", "CC", "CG", "CT", "TA", "TG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    //     {'M', {""}},                                       // {"ATG"}
    //     {'N', {"C", "T"}},                                 // {"AAT", "AAC"}
    //     {'P', {"A", "C", "G", "T"}},                        // {"CCT", "CCC", "CCA", "CCG"}
    //     {'Q', {"A", "G"}},                                 // {"CAA", "CAG"}
    //     {'R', {"CT", "CC", "CA", "CG", "AA", "AG"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
    //     {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
    //     {'T', {"A", "C", "G", "T"}},                     // {"ACT", "ACC", "ACA", "ACG"}
    //     {'V', {"A", "C", "G", "T"}},                     // {"GTT", "GTC", "GTA", "GTG"}},
    //     {'W', {""}},                                       // {"TGG"}
    //     {'Y', {"C", "T"}},                                 // {"TAT", "TAC"}},
    //     {'*', {"AA", "AG", "GA"}}                          // {"TAA", "TAG", "TGA"}
    // };

    std::unordered_map<char, std::vector<std::string>> reduced_codon_table = {
        {'A', {"GCA", "GCC", "GCG", "GCT"}},                       // {"GCT", "GCC", "GCA", "GCG"}
        {'C', {"TGC", "TGT"}},                                 // {"TGT", "TGC"}
        {'D', {"GAC", "GAT"}},                                 // {"GAT", "GAC"}
        {'E', {"GAA", "GAG"}},                                 // {"GAA", "GAG"}
        {'F', {"TTC", "TTT"}},                                 // {"TTT", "TTC"}
        {'G', {"GGA", "GGC", "GGG", "GGT"}},                       // {"GGT", "GGC", "GGA", "GGG"}
        {'H', {"CAC", "CAT"}},           // {"CAT", "CAC"}
        {'I', {"ATA", "ATC", "ATT"}},                            // {"ATT", "ATC", "ATA"}
        {'K', {"AAA", "AAG"}},                                 // {"AAA", "AAG"}
        {'L', {"CTA", "CTC", "CTG", "CTT", "TTA", "TTG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', {"ATG"}},                                       // {"ATG"}
        {'N', {"AAC", "AAT"}},                                 // {"AAT", "AAC"}
        {'P', {"CCA", "CCC", "CCG", "CCT"}},                       // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', {"CAA", "CAG"}},                                 // {"CAA", "CAG"}
        {'R', {"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', {"ACA", "ACC", "ACG", "ACT"}},                       // {"ACT", "ACC", "ACA", "ACG"}
        {'V', {"GTA", "GTC", "GTG", "GTT"}},                       // {"GTT", "GTC", "GTA", "GTG"}},
        {'W', {"TGG"}},                                       // {"TGG"}
        {'Y', {"TAC", "TAT"}},                                 // {"TAT", "TAC"}},
        {'*', {"TAA", "TAG", "TGA"}}                          // {"TAA", "TAG", "TGA"}
    };

    // std::unordered_map<char, int> reduced_codon_lengths_table = {
    //     {'A', 1},
    //     {'C', 1},
    //     {'D', 1},
    //     {'E', 1},
    //     {'F', 1},
    //     {'G', 1},
    //     {'H', 1},
    //     {'I', 1},
    //     {'K', 1},
    //     {'L', 2},
    //     {'M', 0},
    //     {'N', 1},
    //     {'P', 1},
    //     {'Q', 1},
    //     {'R', 2},
    //     {'S', 3},
    //     {'T', 1},
    //     {'V', 1},
    //     {'W', 0},
    //     {'Y', 1},
    //     {'*', 2}};

    std::unordered_map<char, int> reduced_codon_lengths_table = {
        {'A', 3},
        {'C', 3},
        {'D', 3},
        {'E', 3},
        {'F', 3},
        {'G', 3},
        {'H', 3},
        {'I', 3},
        {'K', 3},
        {'L', 3},
        {'M', 3},
        {'N', 3},
        {'P', 3},
        {'Q', 3},
        {'R', 3},
        {'S', 3},
        {'T', 3},
        {'V', 3},
        {'W', 3},
        {'Y', 3},
        {'*', 3}};

    std::unordered_map<std::string, char> translate_codon_table = {
        // Phenylalanine
        {"TTT", 'F'},
        {"TTC", 'F'},
        // Leucine
        {"TTA", 'L'},
        {"TTG", 'L'},
        {"CTT", 'L'},
        {"CTC", 'L'},
        {"CTA", 'L'},
        {"CTG", 'L'},
        // Isoleucine
        {"ATT", 'I'},
        {"ATC", 'I'},
        {"ATA", 'I'},
        // Methionine (Start)
        {"ATG", 'M'},
        // Valine
        {"GTT", 'V'},
        {"GTC", 'V'},
        {"GTA", 'V'},
        {"GTG", 'V'},
        // Serine
        {"TCT", 'S'},
        {"TCC", 'S'},
        {"TCA", 'S'},
        {"TCG", 'S'},
        {"AGT", 'S'},
        {"AGC", 'S'},
        // Proline
        {"CCT", 'P'},
        {"CCC", 'P'},
        {"CCA", 'P'},
        {"CCG", 'P'},
        // Threonine
        {"ACT", 'T'},
        {"ACC", 'T'},
        {"ACA", 'T'},
        {"ACG", 'T'},
        // Alanine
        {"GCT", 'A'},
        {"GCC", 'A'},
        {"GCA", 'A'},
        {"GCG", 'A'},
        // Tyrosine
        {"TAT", 'Y'},
        {"TAC", 'Y'},
        // Histidine
        {"CAT", 'H'},
        {"CAC", 'H'},
        // Glutamine
        {"CAA", 'Q'},
        {"CAG", 'Q'},
        // Asparagine
        {"AAT", 'N'},
        {"AAC", 'N'},
        // Lysine
        {"AAA", 'K'},
        {"AAG", 'K'},
        // Aspartic Acid
        {"GAT", 'D'},
        {"GAC", 'D'},
        // Glutamic Acid
        {"GAA", 'E'},
        {"GAG", 'E'},
        // Cysteine
        {"TGT", 'C'},
        {"TGC", 'C'},
        // Tryptophan
        {"TGG", 'W'},
        // Arginine
        {"CGT", 'R'},
        {"CGC", 'R'},
        {"CGA", 'R'},
        {"CGG", 'R'},
        {"AGA", 'R'},
        {"AGG", 'R'},
        // Glycine
        {"GGT", 'G'},
        {"GGC", 'G'},
        {"GGA", 'G'},
        {"GGG", 'G'},
        // Stop codons
        {"TAA", '*'},
        {"TAG", '*'},
        {"TGA", '*'}};
};

class SIRIUSInstance
{
public:
    int n;
    int max_priority;
    int decidable_protein_length;
    int warm_start_largest_fragment_length;
    SIRIUSTables tables;
    std::string dna_with_holes;
    std::string target_protein;
    std::string decidable_protein;
    std::vector<std::vector<std::string>> warm_start_solution_seqs;

    SIRIUSInstance(
        int n,
        std::string target_protein,
        SIRIUSTables tables,
        const std::string& warm_start_path = "")
        : n(n), target_protein(std::move(target_protein)), tables(std::move(tables))
    {
        if (!warm_start_path.empty()) {
            this->warm_start_solution_seqs = process_warm_start_file(
                warm_start_path,
                this->tables.invariant_codon_table,
                this->tables.translate_codon_table,
                this->tables.skip_aa);

            std::vector<std::vector<std::string>> reduced_warm_start_solution_seqs = process_warm_start_file(
                warm_start_path,
                this->tables.invariant_codon_table,
                this->tables.translate_codon_table,
                this->tables.skip_aa);

            this->warm_start_largest_fragment_length = find_max_fragment_length_from_seqs(reduced_warm_start_solution_seqs);
        }

        this->max_priority = 0;
        this->dna_with_holes = "";
        this->decidable_protein = "";
        this->decidable_protein_length = 0;

        for (const char c : this->target_protein)
        {
            max_priority += this->tables.reduced_codon_lengths_table[c];
            dna_with_holes += this->tables.invariant_codon_table[c];

            if (this->tables.skip_aa.find(c) == this->tables.skip_aa.end())
            {   // if c not in skip_aa
                decidable_protein += c;
                ++decidable_protein_length;
            }

            // decidable_protein += c;
            // ++decidable_protein_length;
        }
    }

    std::vector<std::vector<std::string>> process_warm_start_file(
        const std::string& filename,
        const std::unordered_map<char, std::string>& codon_table,
        const std::unordered_map<std::string, char>& translate_codon_table,
        const std::set<char>& skip_aa_table,
        int codon_length = 3)
    {
        std::string line;
        std::ifstream infile(filename);
        std::vector<std::vector<std::string>> seqs;
    
        // Read and process each warm start sequence
        while (std::getline(infile, line)) {
            std::string protein_translation = translate_dna_to_protein(line, tables);
            if (protein_translation != target_protein) {
                std::cerr << "Error: sequence '" << line
                          << "' translates to '" << protein_translation
                          << "' which does not match target protein '" << target_protein << "'\n";
                std::exit(EXIT_FAILURE);
            }

            std::vector<std::string> variable_bases_of_codons;
            
            for (int i = 0; i + (codon_length - 1) < line.size(); i += codon_length)
            {
                std::string codon = line.substr(i, codon_length);
                char aa = translate_codon_table.at(codon);
                
                if (skip_aa_table.find(aa) != skip_aa_table.end())
                {
                    continue;
                }
                
                std::string bases = codon_table.at(aa);
                std::string variable_bases = "";

                for (int base_i = 0; base_i < bases.size(); ++base_i)
                {
                    if (bases.at(base_i) == '_')
                    {
                        variable_bases += codon.at(base_i);
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

    int find_max_fragment_length_from_seqs(std::vector<std::vector<std::string>> seqs)
    {
        int max_length = 0;

        // Flatten each warm-start sequence (vector<string> --> string)
        std::vector<std::string> flat_seqs;
        for (const auto& seq : seqs) {
            std::string flat;
            for (const std::string& token : seq) {
                flat += token;
            }
            flat_seqs.push_back(flat);
        }

        // std::cout << flat_seqs[0] << std::endl;

        // Compare every unique pair of sequences
        for (size_t i = 0; i < flat_seqs.size(); ++i) {
            for (size_t j = i + 1; j < flat_seqs.size(); ++j) {
                const std::string& s1 = flat_seqs[i];
                const std::string& s2 = flat_seqs[j];

                int current = 0;
                int local_max = 0;
                for (size_t k = 0; k < std::min(s1.size(), s2.size()); ++k) {
                    if (s1[k] == s2[k]) {
                        ++current;
                        local_max = std::max(local_max, current);
                    } else {
                        current = 0;
                    }
                }

                max_length = std::max(max_length, local_max);
            }
        }

        return max_length;
    }

};

// https://github.com/google/or-tools/blob/stable/ortools/sat/sat_parameters.proto
class SIRIUSConfig
{
public:
    operations_research::sat::SatParameters parameters;

    SIRIUSConfig(
        bool log_search_progress,
        unsigned int num_workers,
        unsigned int linearization_level,
        double relative_gap_limit,
        bool fix_variables_to_their_hinted_value = false)
    {
        parameters.set_num_workers(num_workers);
        parameters.set_log_search_progress(log_search_progress);
        parameters.set_linearization_level(linearization_level);
        parameters.set_relative_gap_limit(relative_gap_limit);
        parameters.set_fix_variables_to_their_hinted_value(fix_variables_to_their_hinted_value);
    }
};

class SIRIUSSolver
{
    public:
    SIRIUSTables tables;
    SIRIUSConfig config;
    SIRIUSInstance instance;
    std::unique_ptr<operations_research::sat::Model> model;
    operations_research::sat::CpModelBuilder cp_model;
    operations_research::sat::CpSolverResponse response;

    int dna_size;
    int codon_len;
    int max_priority;
    int target_protein_length;
    int decidable_protein_length;
    std::string dna_with_holes;
    std::string decidable_protein;
    std::unordered_map<std::string, int> map_var_name_to_val;
    VectorBoolVariables all_vars;
    // VectorBoolVariables all_mults;
    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;
    VectorBoolVariables obj_mults;
    std::vector<operations_research::sat::IntVar> additive_obj_mults;
    std::vector<int> prev_obj_vals;
    VectorBoolVariables constraint_obj_mults;
    VectorColumnYVars all_pairs_y_terms;
    VectorColumnZVars all_pairs_z_terms;
    VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;
    std::vector<std::vector<VectorBoolVariables>> all_pairs_all_priorities_combos;
    operations_research::sat::LinearExpr objective;

    SIRIUSSolver(SIRIUSInstance instance, SIRIUSConfig config, SIRIUSTables tables)
        : instance(std::move(instance)),
          config(std::move(config)),
          tables(std::move(tables)),
          model(std::make_unique<operations_research::sat::Model>()),
          cp_model()
    {
        this->codon_len = 3;
        this->max_priority = 0;
        this->decidable_protein_length = 0;
        this->dna_with_holes = "";

        this->dna_size = this->instance.target_protein.size() * 3;
        this->target_protein_length = this->instance.target_protein.size();
        this->max_priority = this->instance.max_priority;
        this->dna_with_holes = this->instance.dna_with_holes;
        this->decidable_protein = this->instance.decidable_protein;
        this->decidable_protein_length = this->instance.decidable_protein_length;
    }

    void assign_var_values_from_solution()
    {
        std::cout << "Assigning var vals from solution...\n";

        int colcount = 0;
        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                for (int aa_pos_i = 0; aa_pos_i < this->decidable_protein_length; ++aa_pos_i)
                {
                    std::string x = this->instance.warm_start_solution_seqs.at(s).at(aa_pos_i);
                    std::string y = this->instance.warm_start_solution_seqs.at(t).at(aa_pos_i);

                    std::string codon_x_bases = "";
                    std::string codon_y_bases = "";

                    std::vector<std::string> codons_vec = this->tables.reduced_codon_table.at(this->decidable_protein.at(aa_pos_i));

                    auto it_x = std::find(codons_vec.begin(), codons_vec.end(), x);
                    auto it_y = std::find(codons_vec.begin(), codons_vec.end(), y);

                    // Calculate the index by subtracting the beginning iterator
                    int codon_x_i = std::distance(codons_vec.begin(), it_x);
                    int codon_y_i = std::distance(codons_vec.begin(), it_y);

                    for (int base_idx = 0; base_idx < this->tables.reduced_codon_lengths_table.at(this->decidable_protein.at(aa_pos_i)); ++base_idx)
                    {
                        std::string var_name_x = absl::StrFormat(
                            "%c%d%d%d%d",
                            codons_vec.at(codon_x_i).at(base_idx),
                            s,
                            aa_pos_i,
                            codon_x_i,
                            base_idx);

                        std::string var_name_y = absl::StrFormat(
                            "%c%d%d%d%d",
                            codons_vec.at(codon_y_i).at(base_idx),
                            t,
                            aa_pos_i,
                            codon_y_i,
                            base_idx);
                        
                        this->map_var_name_to_val.at(var_name_x) = 1;
                        this->map_var_name_to_val.at(var_name_y) = 1;

                        codon_x_bases += var_name_x;
                        codon_y_bases += var_name_y;

                        if (codons_vec.at(codon_x_i).at(base_idx) == codons_vec.at(codon_y_i).at(base_idx))
                        {
                            // y - bind the same bases in different sequences together if they're the same
                            std::string mult_var_name = var_name_x + var_name_y;
                            this->map_var_name_to_val.at(mult_var_name) = 1;

                            // z - if on a column, any y == 1, then the z for that column is 1 as well - easy.
                            std::string z_var_name = absl::StrFormat("z%d", colcount);
                            this->map_var_name_to_val.at(z_var_name) = 1;
                        }

                        ++colcount;
                    }

                    // bind a codon's bases together
                    this->map_var_name_to_val.at(codon_x_bases) = 1;
                    this->map_var_name_to_val.at(codon_y_bases) = 1;
                }
            }
        }
    }

    void set_larger_fragment_obj_vals_to_zero()
    {
        std::cout << "Set larger fragment obj vals to zero...\n";
        for (int i = this->instance.warm_start_largest_fragment_length + 1; i < this->max_priority; ++i)
        {
            this->prev_obj_vals.push_back(0);  // for invariant bases, there will be a default obj value 
                                                // we cannot go under. 
        }
    }

    void set_all_var_vals_to_zero()
    {
        std::cout << "Setting all var vals to zero...\n";
        for (const operations_research::sat::BoolVar &v : this->all_vars)
        {
            this->map_var_name_to_val[v.Name()] = 0;
        }
    }

    void init_new_model()
    {
        std::cout << "Init new model...\n";
        all_vars.clear();

        // TODO model can be cloned
        // https://github.com/google/or-tools/blob/stable/ortools/sat/docs/model.md#introduction
        cp_model = operations_research::sat::CpModelBuilder();
        model = std::make_unique<operations_research::sat::Model>();
        model->Add(operations_research::sat::NewSatParameters(this->config.parameters));
    }

    void build_model(const int current_priority)
    {
        std::cout << "Building model...\n";

        std::cout << "Creating base vars...\n";
        this->sequence_vars_list = add_base_variables();

        std::cout << "Creating codon constraints...\n";
        this->sequence_codons_list = add_codon_constraints();

        std::cout << "Creating symmetry breaking constraints...\n";
        add_codon_mult_relaxed_anti_symmetry_constraints();

        std::cout << "Creating Y chained vars...\n";
        this->all_pairs_y_terms = create_y_chained_vars();

        std::cout << "Creating Z chained vars...\n";
        this->all_pairs_z_terms = create_z_chained_vars();

        std::cout << "Creating lex objective...\n";
        this->obj_mults = create_lex_objective(current_priority);

        std::cout << "Creating objectives and constraints...\n";
        this->constraint_obj_mults = create_forward_chain_lex_objectives_as_constraints(current_priority);
    }

    void set_minimize_objective_value()
    {
        std::cout << "Setting minimize objective...\n";
        this->objective = operations_research::sat::LinearExpr::Sum(this->obj_mults);
        this->cp_model.Minimize(this->objective);
    }

    void add_hints()
    {
        std::cout << "Adding hints from prev solve...\n";
        for (const operations_research::sat::BoolVar &v : all_vars)
        {
            if (this->map_var_name_to_val.find(v.Name()) != this->map_var_name_to_val.end())
            {
                this->cp_model.AddHint(v, this->map_var_name_to_val.at(v.Name()));
            }
        }
    }

    void solve_model()
    {
        std::cout << "Solving model...\n";
        // Solve the model
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solving...\n";
        this->response = SolveCpModel(cp_model.Build(), model.get());
        check_response(this->response);
    }

    void store_solution()
    {
        std::cout << "Storing solution...\n";
        this->map_var_name_to_val.clear();

        for (const operations_research::sat::BoolVar &v : all_vars)
        {
            this->map_var_name_to_val[v.Name()] = operations_research::sat::SolutionIntegerValue(this->response, v);
        }

        this->prev_obj_vals.push_back(this->response.objective_value());
    }

    VectorSequenceAminoacidCodonBoolVariables add_base_variables()
    {
        // Variables for each base
        VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;

        for (int sequence_n = 0; sequence_n < this->instance.n; ++sequence_n)
        {
            VectorAminoacidCodonBoolVariables this_sequence_vars_list_of_list;

            for (int amino_acid_position = 0; amino_acid_position < this->decidable_protein_length; ++amino_acid_position)
            {
                char amino_acid = this->decidable_protein[amino_acid_position];
                VectorCodonBoolVariables codon_vars_list;

                for (size_t codon_number = 0; codon_number < this->tables.reduced_codon_table.at(amino_acid).size(); ++codon_number)
                {
                    const std::string &codon = this->tables.reduced_codon_table.at(amino_acid).at(codon_number);
                    VectorBoolVariables base_vars_list;

                    for (size_t base_idx = 0; base_idx < this->tables.reduced_codon_lengths_table.at(amino_acid); ++base_idx)
                    {
                        std::string var_name = absl::StrFormat(
                            "%c%d%d%d%d",
                            codon[base_idx],
                            sequence_n,
                            amino_acid_position,
                            codon_number,
                            base_idx);

                        operations_research::sat::BoolVar new_bool_var = this->cp_model.NewBoolVar().WithName(var_name);
                        // std::cout << "Creating var " << var_name << "\n";

                        // this->tables.variable_codon_table.find()

                        this->all_vars.push_back(new_bool_var);
                        base_vars_list.push_back(new_bool_var);
                    }
                    codon_vars_list.push_back(base_vars_list);
                }
                this_sequence_vars_list_of_list.push_back(codon_vars_list);
            }
            sequence_vars_list.push_back(this_sequence_vars_list_of_list);
        }

        return sequence_vars_list;
    }

    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> add_codon_constraints()
    {
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;

        // Constraints for Valid Codons
        for (const auto &this_sequence_vars_list_of_list_it : this->sequence_vars_list)
        {
            std::vector<std::vector<operations_research::sat::BoolVar>> this_sequence_codons;

            for (const auto &codon_vars_list_of_lists : this_sequence_vars_list_of_list_it)
            {
                VectorBoolVariables codon_mult_list;

                for (const auto &codon_vars_list : codon_vars_list_of_lists)
                {
                    std::string var_name;
                    for (const auto &var : codon_vars_list)
                    {
                        var_name += var.Name();
                    }

                    // std::cout << var_name << std::endl;
                    operations_research::sat::BoolVar mult = this->cp_model.NewBoolVar().WithName(var_name);

                    const int group_size = codon_vars_list.size();
                    operations_research::sat::LinearExpr group_sum = operations_research::sat::LinearExpr::Sum(codon_vars_list);
                    // Enforce bi-directional implication
                    this->cp_model.AddEquality(group_sum, group_size).OnlyEnforceIf(mult);
                    this->cp_model.AddLessThan(group_sum, group_size).OnlyEnforceIf(mult.Not());

                    this->cp_model.AddEquality(operations_research::sat::LinearExpr::Sum(codon_vars_list), codon_vars_list.size() * mult);

                    this->all_vars.push_back(mult);
                    codon_mult_list.push_back(mult);
                }

                if (!codon_mult_list.empty())
                {
                    // Ensure only one codon auxiliary variable is chosen
                    this->cp_model.AddExactlyOne(codon_mult_list);
                    this_sequence_codons.push_back(codon_mult_list);
                }
            }

            sequence_codons_list.push_back(this_sequence_codons);
        }

        return sequence_codons_list;
    }

    void add_codon_mult_anti_symmetry_constraints()
    {
        int n_sequences = this->sequence_codons_list.size();

        for (int k = 0; k < n_sequences - 1; ++k) {
            int seq_length = this->sequence_codons_list[k].size();
            std::vector<operations_research::sat::BoolVar> prefix_equal(seq_length);
    
            for (int pos = 0; pos < seq_length; ++pos) {
                // Determine if seq[k] and seq[k+1] codons at this position are equal
                prefix_equal[pos] = this->cp_model.NewBoolVar();
                int num_codons = this->sequence_codons_list[k][pos].size();
    
                // Equality at this position across all codon variables
                std::vector<operations_research::sat::BoolVar> codon_equals;
                for (int c = 0; c < num_codons; ++c) {
                    operations_research::sat::BoolVar codon_match = cp_model.NewBoolVar();
                    this->cp_model.AddEquality(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(codon_match);
                    this->cp_model.AddNotEqual(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(codon_match.Not());
                    codon_equals.push_back(codon_match);
                }
                // All codon bits must match to declare equality at this position
                this->cp_model.AddBoolAnd(codon_equals).OnlyEnforceIf(prefix_equal[pos]);

                std::vector<operations_research::sat::BoolVar> codon_differs;
                for (const auto& eq : codon_equals) {
                    codon_differs.push_back(eq.Not());
                }
                this->cp_model.AddBoolOr(codon_differs).OnlyEnforceIf(prefix_equal[pos].Not());

                operations_research::sat::BoolVar lex_less_or_equal = this->cp_model.NewBoolVar();

                // Ensures prefix_equal[pos] = AND(codon_equals)
    
                // Enforce lex ordering at this position
                // Create ordering constraint: seq[k] <= seq[k+1]
                this->cp_model.AddImplication(prefix_equal[pos].Not(), lex_less_or_equal);

                bool jumped = false;
                int c_wrapping_indexer = k % num_codons;

                for (int c = c_wrapping_indexer; c < num_codons;)
                {
                    if (!jumped)
                    {
                        if (c + 1 >= num_codons)
                        {
                            this->cp_model.AddGreaterThan(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(lex_less_or_equal);
                            c = 0;
                            this->cp_model.AddLessThan(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(lex_less_or_equal);
                            ++c;
                        }
                        else
                        {
                            for (int up_to_c = 0; up_to_c < c; ++up_to_c)
                            {
                                this->cp_model.AddLessOrEqual(this->sequence_codons_list[k][pos][up_to_c], this->sequence_codons_list[k+1][pos][up_to_c]).OnlyEnforceIf(lex_less_or_equal);
                            }
                            
                            this->cp_model.AddGreaterThan(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(lex_less_or_equal);
                            this->cp_model.AddLessThan(this->sequence_codons_list[k][pos][c+1], this->sequence_codons_list[k+1][pos][c+1]).OnlyEnforceIf(lex_less_or_equal);
                            c += 2;
                        }

                        jumped = true;
                    }
                    else
                    {
                        this->cp_model.AddLessOrEqual(this->sequence_codons_list[k+1][pos][c], this->sequence_codons_list[k][pos][c]).OnlyEnforceIf(lex_less_or_equal);
                        ++c;
                    }
                }
            }
        }
    }

    // Function courtesy of Google Gemini because there was steam coming out of my head figuring this out.
    //
    // Function to generate and print the sequence based on given inputs
    // n_sequences: The total desired length of the output sequence.
    // num_codons: The base length of the repeating pattern (e.g., 4 for '1 2 3 4').
    std::vector<int> generateSequence(int n_sequences, int num_codons) {
        std::vector<int> patterns;

        // Calculate the number of full cycles of '1' to 'num_codons' that fit into n_sequences.
        // For example, if n_sequences = 10 and num_codons = 4, num_full_cycles will be 2 (10 / 4 = 2).
        int num_full_cycles = n_sequences / num_codons;

        // Calculate the number of remaining elements after the full cycles.
        // For example, if n_sequences = 10 and num_codons = 4, remainder_elements will be 2 (10 % 4 = 2).
        int remainder_elements = n_sequences % num_codons;

        // Loop through each position in the output sequence, from 0 up to n_sequences - 1.
        // The loop iterates exactly n_sequences times to produce the desired number of elements.
        for (int i = 0; i < n_sequences; ++i) {
            int value_to_print; // Variable to store the calculated value for the current position

            // Check if the current index 'i' falls within the range of the full cycles.
            // This handles segments like '1 2 3 4' and '1 2 3 4' for n_sequences = 10, num_codons = 4.
            if (i < num_full_cycles * num_codons) {
                // For full cycles, the value is simply the current index modulo num_codons, plus 1.
                // (i % num_codons) gives a 0-indexed value (0 to num_codons-1),
                // adding 1 converts it to the 1-indexed values (1 to num_codons) shown in the output.
                value_to_print = (i % num_codons) + 1;
            } else {
                // This 'else' block is executed only for the remaining elements after the full cycles.
                // This part handles the "tail" pattern, like '3 4' for n_sequences = 10, num_codons = 4.

                // Calculate the current position within this remainder segment (0-indexed).
                // For example, if i=8 and num_full_cycles*num_codons = 8, current_k_in_remainder will be 0.
                int current_k_in_remainder = i - (num_full_cycles * num_codons);

                // Determine the starting offset for the cyclic pattern within the 1-to-num_codons sequence.
                // This offset is derived from the observed output patterns in the image:
                // - If 1 remainder element (e.g., n_sequences=5, num_codons=4), the offset is (4-1)=3,
                //   meaning the cycle starts from the 4th element (value 4).
                // - If 2 remainder elements (e.g., n_sequences=6, num_codons=4), the offset is (4-2)=2,
                //   meaning the cycle starts from the 3rd element (value 3).
                // - If 3 remainder elements (e.g., n_sequences=7, num_codons=4), the offset is (4-3)=1,
                //   meaning the cycle starts from the 2nd element (value 2).
                // This 'start_offset_for_remainder_0_indexed' represents the 0-indexed position
                // within a '0, 1, ..., num_codons-1' sequence where the remainder segment effectively begins.
                int start_offset_for_remainder_0_indexed = (num_codons - remainder_elements);

                // Calculate the final value to print for the current position in the remainder.
                // This uses the 'start_offset_for_remainder_0_indexed' and 'current_k_in_remainder'
                // to cyclically pick elements from the 1-to-num_codons range.
                // The modulo operator ensures the cycle repeats correctly.
                // Adding 1 converts the 0-indexed result to the desired 1-indexed output values.
                value_to_print = (start_offset_for_remainder_0_indexed + current_k_in_remainder) % num_codons + 1;
            }

            // Print the calculated value for the current position.
            // std::cout << value_to_print;
            patterns.push_back(value_to_print);

            // Print a space after the number, unless it's the very last number in the sequence.
            // This ensures numbers are separated by spaces but there isn't a trailing space.
            // if (i < n_sequences - 1) {
            //     std::cout << " ";
            // }
        }
        // After printing all numbers, print a newline character to move to the next line
        // for any subsequent output, ensuring clean formatting.
        // std::cout << std::endl;
        return patterns;
    }

    void add_codon_mult_relaxed_anti_symmetry_constraints()
    {
        int n_sequences = this->sequence_codons_list.size();

        for (int k = 0; k < n_sequences - 1; ++k) {
            int seq_length = this->sequence_codons_list[k].size();
            std::vector<operations_research::sat::BoolVar> prefix_equal(seq_length);
    
            for (int pos = 0; pos < seq_length; ++pos) {
                // Determine if seq[k] and seq[k+1] codons at this position are equal
                prefix_equal[pos] = this->cp_model.NewBoolVar();
                int num_codons = this->sequence_codons_list[k][pos].size();

                std::vector<int> pattern_end_ranges = generateSequence(n_sequences, num_codons);
    
                // Equality at this position across all codon variables
                std::vector<operations_research::sat::BoolVar> codon_equals;
                for (int c = 0; c < num_codons; ++c) {
                    operations_research::sat::BoolVar codon_match = cp_model.NewBoolVar();
                    this->cp_model.AddEquality(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(codon_match);
                    this->cp_model.AddNotEqual(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k+1][pos][c]).OnlyEnforceIf(codon_match.Not());
                    codon_equals.push_back(codon_match);
                }
                // All codon bits must match to declare equality at this position
                this->cp_model.AddBoolAnd(codon_equals).OnlyEnforceIf(prefix_equal[pos]);

                std::vector<operations_research::sat::BoolVar> codon_differs;
                for (const auto& eq : codon_equals) {
                    codon_differs.push_back(eq.Not());
                }
                this->cp_model.AddBoolOr(codon_differs).OnlyEnforceIf(prefix_equal[pos].Not());

                operations_research::sat::BoolVar lex_less_or_equal = this->cp_model.NewBoolVar();

                // Ensures prefix_equal[pos] = AND(codon_equals)
    
                // Enforce lex ordering at this position
                // Create ordering constraint: seq[k] <= seq[k+1]
                this->cp_model.AddImplication(prefix_equal[pos].Not(), lex_less_or_equal);                
                std::vector<operations_research::sat::BoolVar> selector_vars;  // one for each (i,j) pair

                int selector_index = 0;

                int i_wrapping_indexer = k % num_codons;
                for (int i = i_wrapping_indexer; i < pattern_end_ranges.at(k); ++i)
                {
                    int j_wrapping_indexer = (i + 1) % num_codons;
                    for (int j = j_wrapping_indexer; j < pattern_end_ranges.at(k + 1); ++j)
                    {
                        operations_research::sat::BoolVar selector = this->cp_model.NewBoolVar();
                        selector_vars.push_back(selector);

                        // Enforce the constraint group only if this selector is active
                        this->cp_model.AddGreaterThan(this->sequence_codons_list[k][pos][i], this->sequence_codons_list[k + 1][pos][i])
                                .OnlyEnforceIf({lex_less_or_equal, selector});

                        this->cp_model.AddLessThan(this->sequence_codons_list[k][pos][j], this->sequence_codons_list[k + 1][pos][j])
                                .OnlyEnforceIf({lex_less_or_equal, selector});

                        for (int bi = 0; bi < num_codons; ++bi) {
                            if (bi != i && bi != j) {
                                this->cp_model.AddLessOrEqual(this->sequence_codons_list[k][pos][bi], this->sequence_codons_list[k + 1][pos][bi])
                                        .OnlyEnforceIf({lex_less_or_equal, selector});
                            }
                        }
                    }
                }

                // Create a dummy variable to accumulate the sum of selector_vars
                operations_research::sat::LinearExpr selector_sum;
                for (const auto& var : selector_vars) {
                    selector_sum += var;
                }

                // selector_sum == 1 <=> exactly one selector is active
                // Enforce this only if lex_less_or_equal is true
                this->cp_model.AddEquality(selector_sum, 1).OnlyEnforceIf(lex_less_or_equal);
            }
        }
    }

    VectorColumnYVars create_y_chained_vars()
    {
        VectorColumnYVars all_pairs_y_terms;

        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                VectorLinearExprs this_pair_y_terms;

                for (int aa_pos_i = 0; aa_pos_i < this->decidable_protein_length; ++aa_pos_i)
                {
                    for (int codon_i = 0; codon_i < this->tables.reduced_codon_lengths_table.at(this->decidable_protein.at(aa_pos_i)); ++codon_i)
                    {
                        operations_research::sat::LinearExpr y_terms;

                        for (size_t codon_position_i = 0; codon_position_i < this->tables.reduced_codon_table.at(this->decidable_protein.at(aa_pos_i)).size(); ++codon_position_i)
                        {
                            operations_research::sat::BoolVar x = this->sequence_vars_list.at(s).at(aa_pos_i).at(codon_position_i).at(codon_i);

                            for (size_t codon_position_j = 0; codon_position_j < this->tables.reduced_codon_table.at(this->decidable_protein.at(aa_pos_i)).size(); ++codon_position_j)
                            {
                                operations_research::sat::BoolVar y = this->sequence_vars_list.at(t).at(aa_pos_i).at(codon_position_j).at(codon_i);

                                std::string var_name = x.Name() + y.Name();
                                operations_research::sat::BoolVar z = cp_model.NewBoolVar().WithName(var_name);
                                
                                if (x.Name().front() == y.Name().front())
                                {
                                    // std::cout << "Creating var " << var_name << "\n";
                                    cp_model.AddEquality(z, y).OnlyEnforceIf(x);
                                    cp_model.AddEquality(z, 0).OnlyEnforceIf(~x);
                                }
                                else
                                {
                                    cp_model.AddEquality(z, 0);
                                }

                                // Add z to y_terms
                                y_terms += z;
                                all_vars.push_back(z);
                            }
                        }
                        this_pair_y_terms.push_back(y_terms);
                    }
                }
                all_pairs_y_terms.push_back(this_pair_y_terms);
            }
        }
        return all_pairs_y_terms;
    }

    VectorColumnZVars create_z_chained_vars()
    {
        VectorColumnZVars all_pairs_z_terms;
        int colcount = 0;

        for (const std::vector<operations_research::sat::LinearExpr> v : this->all_pairs_y_terms)
        {
            VectorBoolVariables this_pair_z_terms;

            for (const operations_research::sat::LinearExpr exp : v)
            {
                std::string var_name = absl::StrFormat("z%d", colcount);
                operations_research::sat::BoolVar z_j = this->cp_model.NewBoolVar().WithName(var_name);
                // std::cout << "Creating var " << var_name << "\n";

                this->cp_model.AddEquality(z_j, exp);

                this->all_vars.push_back(z_j);
                this_pair_z_terms.push_back(z_j);
                colcount += 1;
            }

            all_pairs_z_terms.push_back(this_pair_z_terms);
        }

        return all_pairs_z_terms;
    }

    VectorBoolVariables create_lex_objective(const int current_priority)
    {
        VectorBoolVariables all_mults;

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
        {
            std::vector<std::vector<operations_research::sat::BoolVar>> combos = generate_combinations(
                this->all_pairs_z_terms.at(i), current_priority);

            for (const std::vector<operations_research::sat::BoolVar> &vi : combos)
            {
                std::string var_name = "";

                for (const operations_research::sat::BoolVar &i : vi)
                {
                    var_name += i.Name();
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                // std::cout << "Creating var " << var_name << "\n";

                // Binds mult to multiple of z's: mult = z1*z2*z3*...
                for (const operations_research::sat::BoolVar &i : vi)
                {
                    cp_model.AddLessOrEqual(mult, i);
                }
                
                cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));
                // ---

                this->all_vars.push_back(mult);
                all_mults.push_back(mult);
            }
        }

        return all_mults;
    }

    VectorBoolVariables create_forward_chain_lex_objectives_as_constraints(const int current_priority)
    {
        int counter = 0;
        VectorBoolVariables constraint_obj_mults;

        for (int go_down = this->max_priority; go_down > current_priority; --go_down)
        {
            VectorBoolVariables obj_mults;

            for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
            {
                // std::cout << "i: " << i << " godown:" << go_down << std::endl;
                std::vector<VectorBoolVariables> combos = generate_combinations(
                    this->all_pairs_z_terms.at(i), go_down);
                // std::cout << "done\n";
                
                for (VectorBoolVariables vi : combos)
                {
                    std::string var_name = "";

                    for (const operations_research::sat::BoolVar &i : vi)
                    {
                        var_name += i.Name();
                    }

                    operations_research::sat::BoolVar mult = this->cp_model.NewBoolVar().WithName(var_name);
                    // std::cout << "Creating var " << var_name << "\n";

                    for (const operations_research::sat::BoolVar &i : vi)
                    {
                        this->cp_model.AddLessOrEqual(mult, i);
                    }
                    this->cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));

                    this->all_vars.push_back(mult);
                    obj_mults.push_back(mult);
                    constraint_obj_mults.push_back(mult);
                }
            }

            if (this->prev_obj_vals.size() > counter)
            {
                // std::cout << "counter at " << counter << std::endl;
                this->cp_model.AddLessOrEqual(operations_research::sat::LinearExpr::Sum(obj_mults), this->prev_obj_vals.at(counter));
                counter++;
            }
        }

        return constraint_obj_mults;
    }

};

std::string translate_dna_to_protein(const std::string &dna, SIRIUSTables &tables);
std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);
std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);
std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);
void print_length_counts(const std::unordered_map<int, int> &length_counts);
void check_response(const operations_research::sat::CpSolverResponse &response);

int main(int argc, char *argv[])
{
    /*
       I do feel like NO ONE CAN SAVE ME I AM SO ALONE
       And yet I CRIED, I called for HELP, FORSAKEN
       But now I know the only way is UNDERSTAND THE LIVING
       OBEY THE RULE OF LIGHT and FACE THE FEAR inside out.
       LOST, I found there a stone erected in line
       WITH ONE OF THE BRIGHTEST STARS OF ALL THE NIGHT SKY VAULT
       And I took my time, took off the moss
       Washed away the dust and gave A NEW LEASE OF LIFE
       Its MYSTICAL FORCE
       I grab it now and praise this lord of earth and stone
       Make passage for souls awaken
       So it returns to WHERE IT'S ALWAYS BEEN with the gods
    */

    // Systematische Identifikation Redundanter, Identisch Uebersetzter Sequenzen
    ::google::InitGoogleLogging("SIRIUS");

    int num_sequences;
    std::string init_target_protein;

    bool use_warmstart;
    std::string warm_start_path = "";

    std::cout << BLUE << "> " << RESET << "Welcome to " << BLUE << "SIRIUS" << RESET << std::endl;

    if (argc < 3)
    {
        std::cout << BLUE << "> " << RESET << "Enter target protein: ";
        std::cin >> init_target_protein;
    
        std::cout << BLUE << "> " << RESET << "Enter number of sequences: ";
        std::cin >> num_sequences;
    }
    else
    {
        init_target_protein = argv[1];
        num_sequences = std::stoi(argv[2]);
    
        for (int i = 3; i < argc; ++i) {
            std::string arg = argv[i];
            std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
            if (arg == "--warmstart=true") {
                use_warmstart = true;
            }
        }
    }

    std::cout << BLUE << "> " << RESET << "Creating " << num_sequences << "x " << init_target_protein << std::endl;

    if (use_warmstart)
    {
        warm_start_path = "warm_start.txt";
        std::cout << BLUE << "> " << RESET << "Using warm start from " << warm_start_path << std::endl;
    }
    
    SIRIUSTables tables;
    SIRIUSConfig config(true, 16, 0, 0);
    SIRIUSInstance instance(num_sequences, init_target_protein, tables, warm_start_path);
    SIRIUSSolver sirius_solver(instance, config, tables);

    int current_priority = instance.warm_start_largest_fragment_length;

    if (use_warmstart)
    {
        sirius_solver.max_priority = current_priority;
    }
    else
    {
        current_priority = sirius_solver.max_priority;
    }
    
    std::cout <<  BLUE << "> " << RESET << "Max fragment length: " << sirius_solver.max_priority << "\n";
    std::cout <<  BLUE << "> " << RESET << "Current fragment length:  " << current_priority << "\n";

    std::queue<char> int_vars;

    if (use_warmstart)
    {
        sirius_solver.init_new_model();
        sirius_solver.set_larger_fragment_obj_vals_to_zero();
        sirius_solver.build_model(current_priority);
        sirius_solver.set_minimize_objective_value();
        sirius_solver.set_all_var_vals_to_zero();
        sirius_solver.assign_var_values_from_solution();
        sirius_solver.solve_model();
        sirius_solver.store_solution();

        --current_priority;
    }
    
    for (current_priority; current_priority > 0; --current_priority)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Setting up fragments of length " << current_priority << "...\n";

        sirius_solver.init_new_model();
        sirius_solver.build_model(current_priority);
        sirius_solver.set_minimize_objective_value();
        sirius_solver.add_hints();
        sirius_solver.solve_model();
        sirius_solver.store_solution();

        // Variables for each base
        // Flatten sequence_vars_list
        std::vector<operations_research::sat::BoolVar> base_vars;
        for (const auto &a : sirius_solver.sequence_vars_list)
            for (const auto &b : a)
                for (const auto &c : b)
                    for (const auto &var : c)
                        base_vars.push_back(var);

        // Clear
        std::queue<char>().swap(int_vars);
        for (const auto &var : base_vars)
        {
            if (operations_research::sat::SolutionIntegerValue(sirius_solver.response, var))
            {
                int_vars.push(var.Name().front());
            }
        }
    }

    // -----

    // Variables for each base
    // Flatten sequence_vars_list
    std::vector<operations_research::sat::BoolVar> base_vars;
    for (const auto &a : sirius_solver.sequence_vars_list)
        for (const auto &b : a)
            for (const auto &c : b)
                for (const auto &var : c)
                    base_vars.push_back(var);

    // Clear
    std::queue<char>().swap(int_vars);
    for (const auto &var : base_vars)
    {
        if (operations_research::sat::SolutionIntegerValue(sirius_solver.response, var))
        {
            int_vars.push(var.Name().front());
        }
    }

    int counter = 0;
    std::string seq;

    std::vector<std::string> all_out_seqs;
    for (int seq_n = 0; seq_n < num_sequences; ++seq_n)
    {
        std::string this_seq = "";
        for (int i = 0; i < sirius_solver.dna_size; ++i)
        {
            if (sirius_solver.dna_with_holes[i] == '_')
            {
                this_seq += int_vars.front();
                int_vars.pop();
            }
            else
            {
                this_seq += sirius_solver.dna_with_holes[i];
            }
        }
        all_out_seqs.push_back(this_seq);
    }

    // File for writing solutions
    std::ofstream out_file("out.txt");
    if (!out_file)
    {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }
    for (auto seq : all_out_seqs)
    {
        std::cout << seq << std::endl;
        out_file << seq << "\n";
    }
    out_file.close();

    for (size_t i = 0; i < all_out_seqs.size(); ++i)
    {
        std::string protein = translate_dna_to_protein(all_out_seqs[i], tables);

        std::cout << "Protein " << i + 1 << ": " << protein << "\n";
        if (protein != init_target_protein)
        {
            std::cout << "SEQ " << i << " " << all_out_seqs[i] << " NOT THE SAME\n";
        }
    }

    auto [all_stretches, length_counts] = find_all_homologous_stretches_and_count_lengths(all_out_seqs);
    print_length_counts(length_counts);

    std::cout << BLUE << "> " RESET << "[" << elapsed_since_start() << "] Done.\n";

    return 0;
}

// ===============================================

// Create sliding window combinations out of vars depending on combination size
// E.g., Input: vars = [z1, z2, z3], combination_size = 2
//       Output: [[z1, z2], [z2, z3]]
std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size)
{
    std::vector<std::vector<operations_research::sat::BoolVar>> combinations;

    if (combination_size <= 0 || combination_size > vars.size())
        return combinations;

    for (size_t i = 0; i <= vars.size() - combination_size; ++i)
    {
        combinations.emplace_back(vars.begin() + i, vars.begin() + i + combination_size);
    }

    return combinations;
}

std::string elapsed_since_start()
{
    static const auto start = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();

    int hours = elapsed / 3600;
    int minutes = (elapsed % 3600) / 60;
    int seconds = elapsed % 60;

    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << hours << ":"
        << std::setw(2) << std::setfill('0') << minutes << ":"
        << std::setw(2) << std::setfill('0') << seconds;

    return oss.str();
}

std::string translate_dna_to_protein(const std::string &dna, SIRIUSTables &tables)
{
    std::string protein;
    for (size_t i = 0; i + 2 < dna.size(); i += 3)
    {
        std::string codon = dna.substr(i, 3);
        auto it = tables.translate_codon_table.find(codon);
        if (it != tables.translate_codon_table.end())
        {
            protein += it->second;
        }
        else
        {
            protein += 'X'; // Unknown codon
        }
    }
    return protein;
}

// -------------------------------------------------
// For output print
// Find stretches of homology between two sequences
std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2)
{
    std::vector<Stretch> stretches;
    int start = -1;

    for (size_t i = 0; i < seq1.size(); ++i)
    {
        if (seq1[i] == seq2[i])
        {
            if (start == -1)
                start = i;
        }
        else
        {
            if (start != -1)
            {
                stretches.emplace_back(start, i - 1, static_cast<int>(i - start));
                start = -1;
            }
        }
    }

    // Check if a stretch ended at the last character
    if (start != -1)
    {
        stretches.emplace_back(start, static_cast<int>(seq1.size() - 1), static_cast<int>(seq1.size() - start));
    }

    // Sort stretches by descending length
    std::sort(stretches.begin(), stretches.end(),
              [](const Stretch &a, const Stretch &b)
              {
                  return std::get<2>(a) > std::get<2>(b);
              });

    return stretches;
}

// Find all homologous stretches across all sequence pairs and count lengths
std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences)
{
    std::map<std::string, std::vector<Stretch>> all_stretches;
    std::unordered_map<int, int> length_counts;

    int pair_index = 1;
    for (size_t i = 0; i < sequences.size(); ++i)
    {
        for (size_t j = i + 1; j < sequences.size(); ++j)
        {
            std::string key = "Pair " + std::to_string(pair_index) +
                              " (Seq " + std::to_string(i + 1) +
                              " vs Seq " + std::to_string(j + 1) + ")";
            std::vector<Stretch> stretches = find_homologous_stretches(sequences[i], sequences[j]);
            all_stretches[key] = stretches;

            for (const auto &stretch : stretches)
            {
                int length = std::get<2>(stretch);
                ++length_counts[length];
            }

            ++pair_index;
        }
    }

    return {all_stretches, length_counts};
}

void print_length_counts(const std::unordered_map<int, int> &length_counts)
{
    // Convert unordered_map to a vector of pairs for sorting
    std::vector<std::pair<int, int>> sorted_counts(length_counts.begin(), length_counts.end());

    // Sort by length in descending order
    std::sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto &a, const auto &b)
              {
                  return a.first > b.first;
              });

    std::cout << "\nLength Counts:\n";
    for (const auto &[length, count] : sorted_counts)
    {
        std::cout << "Length " << length << ": " << count << " occurrences\n";
    }
}
// -------------------------------------------------

// Print out solver response message
void check_response(const operations_research::sat::CpSolverResponse &response)
{
    if (response.status() == operations_research::sat::CpSolverStatus::OPTIMAL)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Optimal solution found.\n";
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Objective value: " << response.objective_value() << std::endl;
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::FEASIBLE)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Feasible solution found.\n";
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Objective value: " << response.objective_value() << std::endl;
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE)
    {
        std::cout << RED << "> " << RESET << "[" << elapsed_since_start() << "] Infeasible.\n";
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID)
    {
        std::cout << RED << "> " << RESET << "[" << elapsed_since_start() << "] MODEL_INVALID.\n";
    }
}
