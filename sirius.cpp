
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
using VectorAminoacidCodonBoolVariables = std::vector<std::vector<VectorBoolVariables>>;
using VectorSequenceAminoacidCodonBoolVariables = std::vector<VectorAminoacidCodonBoolVariables>;
class SIRIUSTables;
class SIRIUSInstance;
class SIRIUSConfig;
class SIRIUSSolver; // Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;
std::string translate_dna_to_protein(const std::string &dna, SIRIUSTables &tables);
std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);
std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);
std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);
void print_length_counts(const std::unordered_map<int, int> &length_counts);
void check_response(const operations_research::sat::CpSolverResponse &response);

struct SIRIUSTables
{
    std::unordered_map<char, std::string> invariant_codon_table = {
        {'A', "GC_"}, // {"GCT", "GCC", "GCA", "GCG"}
        {'C', "TG_"}, // {"TGT", "TGC"}
        {'D', "GA_"}, // {"GAT", "GAC"}
        {'E', "GA_"}, // {"GAA", "GAG"}
        {'F', "TT_"}, // {"TTT", "TTC"}
        {'G', "GG_"}, // {"GGT", "GGC", "GGA", "GGG"}
        {'H', "CA_"}, // {"CAT", "CAC"}
        {'I', "AT_"}, // {"ATT", "ATC", "ATA"}
        {'K', "AA_"}, // {"AAA", "AAG"}
        {'L', "_T_"}, // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', "ATG"},
        {'N', "AA_"}, // {"AAT", "AAC"}
        {'P', "CC_"}, // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', "CA_"}, // {"CAA", "CAG"}
        {'R', "_G_"}, // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', "___"}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', "AC_"}, // {"ACT", "ACC", "ACA", "ACG"}
        {'V', "GT_"}, // {"GTT", "GTC", "GTA", "GTG"}
        {'W', "TGG"}, // TGG
        {'Y', "TA_"}, // {"TAT", "TAC"}
        {'*', "T__"}  // {"TAA", "TAG", "TGA"}
    };

    // std::unordered_map<char, std::string> invariant_codon_table = {
    //     {'A', "___"}, // {"GCT", "GCC", "GCA", "GCG"}
    //     {'C', "___"}, // {"TGT", "TGC"}
    //     {'D', "___"}, // {"GAT", "GAC"}
    //     {'E', "___"}, // {"GAA", "GAG"}
    //     {'F', "___"}, // {"TTT", "TTC"}
    //     {'G', "___"}, // {"GGT", "GGC", "GGA", "GGG"}
    //     {'H', "___"}, // {"CAT", "CAC"}
    //     {'I', "___"}, // {"ATT", "ATC", "ATA"}
    //     {'K', "___"}, // {"AAA", "AAG"}
    //     {'L', "___"}, // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    //     {'M', "___"},
    //     {'N', "___"}, // {"AAT", "AAC"}
    //     {'P', "___"}, // {"CCT", "CCC", "CCA", "CCG"}
    //     {'Q', "___"}, // {"CAA", "CAG"}
    //     {'R', "___"}, // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
    //     {'S', "___"}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
    //     {'T', "___"}, // {"ACT", "ACC", "ACA", "ACG"}
    //     {'V', "___"}, // {"GTT", "GTC", "GTA", "GTG"}
    //     {'W', "___"}, // TGG
    //     {'Y', "___"}, // {"TAT", "TAC"}
    //     {'*', "___"}  // {"TAA", "TAG", "TGA"}
    // };

    std::set<char> skip_aa = {'M', 'W'};
    // std::set<char> skip_aa = {};

    std::unordered_map<char, std::vector<std::string>> reduced_codon_table = {
        {'A', {"A", "C", "G", "T"}},                       // {"GCT", "GCC", "GCA", "GCG"}
        {'C', {"C", "T"}},                                 // {"TGT", "TGC"}
        {'D', {"C", "T"}},                                 // {"GAT", "GAC"}
        {'E', {"A", "G"}},                                 // {"GAA", "GAG"}
        {'F', {"C", "T"}},                                 // {"TTT", "TTC"}
        {'G', {"A", "C", "G", "T"}},                       // {"GGT", "GGC", "GGA", "GGG"}
        {'H', {"C", "T"}},                                 // {"CAT", "CAC"}
        {'I', {"A", "C", "T"}},                            // {"ATT", "ATC", "ATA"}
        {'K', {"A", "G"}},                                 // {"AAA", "AAG"}
        {'L', {"CA", "CC", "CG", "CT", "TA", "TG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', {""}},                                       // {"ATG"}
        {'N', {"C", "T"}},                                 // {"AAT", "AAC"}
        {'P', {"A", "C", "G", "T"}},                       // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', {"A", "G"}},                                 // {"CAA", "CAG"}
        {'R', {"AA", "AG", "CA", "CC", "CG", "CT"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', {"A", "C", "G", "T"}},                       // {"ACT", "ACC", "ACA", "ACG"}
        {'V', {"A", "C", "G", "T"}},                       // {"GTT", "GTC", "GTA", "GTG"}},
        {'W', {""}},                                       // {"TGG"}
        {'Y', {"C", "T"}},                                 // {"TAT", "TAC"}},
        {'*', {"AA", "AG", "GA"}}                          // {"TAA", "TAG", "TGA"}
    };

    // std::unordered_map<char, std::vector<std::string>> reduced_codon_table = {
    //     {'A', {"GCA", "GCC", "GCG", "GCT"}},                       // {"GCT", "GCC", "GCA", "GCG"}
    //     {'C', {"TGC", "TGT"}},                                 // {"TGT", "TGC"}
    //     {'D', {"GAC", "GAT"}},                                 // {"GAT", "GAC"}
    //     {'E', {"GAA", "GAG"}},                                 // {"GAA", "GAG"}
    //     {'F', {"TTC", "TTT"}},                                 // {"TTT", "TTC"}
    //     {'G', {"GGA", "GGC", "GGG", "GGT"}},                       // {"GGT", "GGC", "GGA", "GGG"}
    //     {'H', {"CAC", "CAT"}},           // {"CAT", "CAC"}
    //     {'I', {"ATA", "ATC", "ATT"}},                            // {"ATT", "ATC", "ATA"}
    //     {'K', {"AAA", "AAG"}},                                 // {"AAA", "AAG"}
    //     {'L', {"CTA", "CTC", "CTG", "CTT", "TTA", "TTG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    //     {'M', {"ATG"}},                                       // {"ATG"}
    //     {'N', {"AAC", "AAT"}},                                 // {"AAT", "AAC"}
    //     {'P', {"CCA", "CCC", "CCG", "CCT"}},                       // {"CCT", "CCC", "CCA", "CCG"}
    //     {'Q', {"CAA", "CAG"}},                                 // {"CAA", "CAG"}
    //     {'R', {"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
    //     {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
    //     {'T', {"ACA", "ACC", "ACG", "ACT"}},                       // {"ACT", "ACC", "ACA", "ACG"}
    //     {'V', {"GTA", "GTC", "GTG", "GTT"}},                       // {"GTT", "GTC", "GTA", "GTG"}},
    //     {'W', {"TGG"}},                                       // {"TGG"}
    //     {'Y', {"TAC", "TAT"}},                                 // {"TAT", "TAC"}},
    //     {'*', {"TAA", "TAG", "TGA"}}                          // {"TAA", "TAG", "TGA"}
    // };

    std::unordered_map<char, int> reduced_codon_lengths_table = {
        {'A', 1},
        {'C', 1},
        {'D', 1},
        {'E', 1},
        {'F', 1},
        {'G', 1},
        {'H', 1},
        {'I', 1},
        {'K', 1},
        {'L', 2},
        {'M', 0},
        {'N', 1},
        {'P', 1},
        {'Q', 1},
        {'R', 2},
        {'S', 3},
        {'T', 1},
        {'V', 1},
        {'W', 0},
        {'Y', 1},
        {'*', 2}};

    // std::unordered_map<char, int> reduced_codon_lengths_table = {
    //     {'A', 3},
    //     {'C', 3},
    //     {'D', 3},
    //     {'E', 3},
    //     {'F', 3},
    //     {'G', 3},
    //     {'H', 3},
    //     {'I', 3},
    //     {'K', 3},
    //     {'L', 3},
    //     {'M', 3},
    //     {'N', 3},
    //     {'P', 3},
    //     {'Q', 3},
    //     {'R', 3},
    //     {'S', 3},
    //     {'T', 3},
    //     {'V', 3},
    //     {'W', 3},
    //     {'Y', 3},
    //     {'*', 3}};

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
        const std::string &warm_start_path = "")
        : n(n),
          target_protein(std::move(target_protein)),
          tables(std::move(tables))
    {
        if (!warm_start_path.empty())
        {
            this->warm_start_solution_seqs = process_warm_start_file(
                warm_start_path,
                this->tables.invariant_codon_table,
                this->tables.translate_codon_table,
                this->tables.skip_aa);

            this->warm_start_largest_fragment_length = find_max_fragment_length_from_seqs(this->warm_start_solution_seqs);
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
            { // if c not in skip_aa
                decidable_protein += c;
                ++decidable_protein_length;
            }
        }
    }

    std::vector<std::vector<std::string>> process_warm_start_file(
        const std::string &filename,
        const std::unordered_map<char, std::string> &codon_table,
        const std::unordered_map<std::string, char> &translate_codon_table,
        const std::set<char> &skip_aa_table,
        int codon_length = 3)
    {
        std::string line;
        std::ifstream infile(filename);
        std::vector<std::vector<std::string>> seqs;

        // Read and process each warm start sequence
        while (std::getline(infile, line))
        {
            std::string protein_translation = translate_dna_to_protein(line, tables);
            if (protein_translation != target_protein)
            {
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

    int find_max_fragment_length_from_seqs(const std::vector<std::vector<std::string>> &seqs)
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
        // parameters.set_max_time_in_seconds(60);
    }
};

class SIRIUSSolver
{
public:
    SIRIUSTables tables;
    SIRIUSConfig config;
    SIRIUSInstance instance;

    operations_research::sat::CpModelBuilder base_cp_model;
    operations_research::sat::CpSolverResponse response;

    int dna_size;
    int max_priority;
    int target_protein_length;
    int decidable_protein_length;
    std::string dna_with_holes;
    std::string decidable_protein;
    VectorBoolVariables all_vars;
    std::vector<std::vector<operations_research::sat::IntVar>> codon_idx;
    std::unordered_map<std::string, int> map_var_name_to_val;

    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;
    std::vector<std::vector<operations_research::sat::BoolVar>> all_obj_mults;
    std::vector<int> all_obj_vals;
    std::vector<int> prev_var_indices;
    std::vector<std::vector<operations_research::sat::BoolVar>> all_pairs_z_terms;
    VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;

    SIRIUSSolver(SIRIUSInstance instance, SIRIUSConfig config, SIRIUSTables tables)
        : instance(std::move(instance)),
          config(std::move(config)),
          tables(std::move(tables))
    {
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
        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                for (int aa_pos_i = 0; aa_pos_i < this->decidable_protein_length; ++aa_pos_i)
                {
                    std::string x = this->instance.warm_start_solution_seqs.at(s).at(aa_pos_i);
                    std::string y = this->instance.warm_start_solution_seqs.at(t).at(aa_pos_i);

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
                    }
                }
            }
        }
    }

    void set_larger_fragment_obj_vals_to_zero()
    {
        std::cout << "Set larger fragment obj vals to zero...\n";
        for (int i = this->instance.warm_start_largest_fragment_length; i < this->all_obj_mults.size() + 1; ++i)
        {
            this->all_obj_vals[i] = 0;
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

    std::vector<int> collect_proto_indices(const VectorBoolVariables &vars)
    {
        std::vector<int> idx;
        for (const operations_research::sat::BoolVar &v : vars)
        {
            idx.push_back(v.index());
        }
        return idx;
    }

    operations_research::sat::CpModelBuilder make_new_cp_model()
    {
        return operations_research::sat::CpModelBuilder();
    }

    operations_research::sat::CpModelBuilder clone_cp_model(operations_research::sat::CpModelBuilder &cp_model)
    {
        return cp_model.Clone();
    }

    void init_cp_model(bool use_warmstart, int priority)
    {
        std::cout << "Building model...\n";
        this->base_cp_model = make_new_cp_model();

        std::cout << "Creating base vars...\n";
        this->sequence_vars_list = add_base_variables(this->base_cp_model);

        std::cout << "Creating codon constraints...\n";
        this->sequence_codons_list = add_codon_constraints(this->base_cp_model);

        std::cout << "Creating Z chained vars...\n";
        this->all_pairs_z_terms = create_z_chained_vars(this->base_cp_model);

        std::cout << "Creating lex objs...\n";
        this->all_obj_vals.resize(priority + 1);
        this->all_obj_mults.resize(priority + 1);
        for (int p = priority; p >= 0; --p)
        {
            this->all_obj_mults[p] = create_lex_objective(this->base_cp_model, p + 1);
        }

        if (use_warmstart)
        {
            set_larger_fragment_obj_vals_to_zero();
        }
    }

    void set_minimize_objective_value(operations_research::sat::CpModelBuilder &cp_model, const int current_priority)
    {
        std::cout << "Setting minimize objective...\n";
        operations_research::sat::LinearExpr obj;
        for (operations_research::sat::BoolVar &v : this->all_obj_mults[current_priority])
        {
            obj += cp_model.GetBoolVarFromProtoIndex(v.index());
        }
        cp_model.Minimize(obj);
    }

   void add_hints(operations_research::sat::CpModelBuilder& cp_model)
    {
        std::cout << "Adding hints from prev solve...\n";
        for (const operations_research::sat::BoolVar& v_old : all_vars)
        {
            auto it = map_var_name_to_val.find(v_old.Name());
            if (it == map_var_name_to_val.end())
            {
                continue;
            }

            // Re-acquire the variable from THIS clone by proto index.
            operations_research::sat::BoolVar v_clone = cp_model.GetBoolVarFromProtoIndex(v_old.index());
            cp_model.AddHint(v_clone, it->second);
        }
    }


    void solve_model(operations_research::sat::CpModelBuilder &cp_model)
    {
        std::cout << "Solving model...\n";
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solving...\n";

        operations_research::sat::Model sat_model;
        sat_model.Add(operations_research::sat::NewSatParameters(config.parameters));
        this->response = SolveCpModel(cp_model.Build(), &sat_model);
        check_response(this->response);
    }

    void store_solution(const int current_priority)
    {
        std::cout << "Storing solution...\n";
        this->map_var_name_to_val.clear();

        for (const operations_research::sat::BoolVar &v : all_vars)
        {
            this->map_var_name_to_val[v.Name()] = operations_research::sat::SolutionIntegerValue(this->response, v);
        }

        this->all_obj_vals[current_priority] = this->response.objective_value();
    }

    std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>>
    add_base_variables(operations_research::sat::CpModelBuilder &cp_model)
    {
        int counter = 0;
        // Variables for each base
        VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;
        this->codon_idx.assign(this->instance.n, std::vector<operations_research::sat::IntVar>(decidable_protein_length));

        for (int sequence_n = 0; sequence_n < this->instance.n; ++sequence_n)
        {
            VectorAminoacidCodonBoolVariables this_sequence_vars_list_of_list;

            for (int amino_acid_position = 0; amino_acid_position < this->decidable_protein_length; ++amino_acid_position)
            {
                char amino_acid = this->decidable_protein[amino_acid_position];
                const std::vector<std::string> &codon_list = this->tables.reduced_codon_table.at(amino_acid);
                int codon_list_size = static_cast<int>(codon_list.size());
                std::vector<std::vector<operations_research::sat::BoolVar>> codon_vars_list;

                for (size_t codon_number = 0; codon_number < codon_list_size; ++codon_number)
                {
                    const std::string &codon = codon_list.at(codon_number);
                    int codon_length = static_cast<int>(codon.length());
                    VectorBoolVariables base_vars_list;

                    for (size_t base_idx = 0; base_idx < codon_length; ++base_idx)
                    {
                        std::string var_name = absl::StrFormat(
                            "%c%d%d%d%d",
                            codon[base_idx],
                            sequence_n,
                            amino_acid_position,
                            codon_number,
                            base_idx);

                        operations_research::sat::BoolVar new_bool_var = cp_model.NewBoolVar().WithName(var_name);
                        ++counter;

                        this->all_vars.push_back(new_bool_var);
                        base_vars_list.push_back(new_bool_var);
                    }
                    codon_vars_list.push_back(base_vars_list);
                }
                this_sequence_vars_list_of_list.push_back(codon_vars_list);
            }
            sequence_vars_list.push_back(this_sequence_vars_list_of_list);
        }

        std::cout << "Created " << counter << " base vars\n";
        return sequence_vars_list;
    }

    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> add_codon_constraints(operations_research::sat::CpModelBuilder &cp_model)
    {
        int counter = 0;
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
                    if (codon_vars_list.size() == 1)
                    {
                        codon_mult_list.push_back(codon_vars_list.at(0));
                        continue;
                    }

                    std::string var_name = "codon_";
                    for (const auto &var : codon_vars_list)
                    {
                        var_name += var.Name();
                    }

                    operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                    ++counter;

                    // mult ⇒ every base var is 1
                    for (const auto& b : codon_vars_list) {
                        cp_model.AddImplication(mult, b);
                    }
                    // (AND bases) ⇒ mult as one clause: (¬b1 ∨ ¬b2 ∨ ... ∨ mult)
                    {
                        std::vector<operations_research::sat::BoolVar> clause;
                        clause.reserve(codon_vars_list.size() + 1);
                        clause.push_back(mult);
                        for (const auto& b : codon_vars_list)
                        {
                            clause.push_back(b.Not());
                        }
                        cp_model.AddBoolOr(clause);
                    }

                    this->all_vars.push_back(mult);
                    codon_mult_list.push_back(mult);
                }

                if (!codon_mult_list.empty())
                {
                    // Ensure only one codon auxiliary variable is chosen
                    cp_model.AddExactlyOne(codon_mult_list);
                    this_sequence_codons.push_back(codon_mult_list);
                }
            }

            sequence_codons_list.push_back(this_sequence_codons);
        }

        std::cout << "Created " << counter << " codon mult vars\n";
        return sequence_codons_list;
    }

    std::vector<std::vector<operations_research::sat::BoolVar>> create_z_chained_vars(operations_research::sat::CpModelBuilder &cp_model)
    {
        int counter = 0;
        std::vector<std::vector<operations_research::sat::BoolVar>> all_pairs_z_terms;

        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                int colcount = 0;
                VectorBoolVariables this_pair_z_terms;

                for (int aa_pos_i = 0; aa_pos_i < this->decidable_protein_length; ++aa_pos_i)
                {
                    for (int codon_i = 0; codon_i < this->tables.reduced_codon_lengths_table.at(this->decidable_protein.at(aa_pos_i)); ++codon_i)
                    {
                        int codon_positions = this->tables.reduced_codon_table.at(this->decidable_protein.at(aa_pos_i)).size();
                        operations_research::sat::LinearExpr z_terms;

                        for (size_t codon_position_i = 0; codon_position_i < codon_positions; ++codon_position_i)
                        {
                            operations_research::sat::BoolVar x = this->sequence_vars_list.at(s).at(aa_pos_i).at(codon_position_i).at(codon_i);

                            for (size_t codon_position_j = 0; codon_position_j < codon_positions; ++codon_position_j)
                            {
                                operations_research::sat::BoolVar y = this->sequence_vars_list.at(t).at(aa_pos_i).at(codon_position_j).at(codon_i);

                                std::string var_name = x.Name() + y.Name();
                                if (x.Name().front() == y.Name().front())
                                {
                                    operations_research::sat::BoolVar z = cp_model.NewBoolVar().WithName(var_name);
                                    ++counter;

                                    cp_model.AddMultiplicationEquality(z, x, y);

                                    // Add z to z_terms
                                    z_terms += z;
                                }
                            }
                        }
                        
                        std::string var_name = absl::StrFormat("z%d_p%d:%d", colcount, s, t);
                        operations_research::sat::BoolVar z_j = cp_model.NewBoolVar().WithName(var_name);
                        colcount += 1;

                        cp_model.AddEquality(z_j, z_terms);
                        this_pair_z_terms.push_back(z_j);
                        // this->all_vars.push_back(z_j);
                    }
                }
                all_pairs_z_terms.push_back(this_pair_z_terms);
            }
        }

        std::cout << "Created " << counter << " y vars\n";
        return all_pairs_z_terms;
    }

    VectorBoolVariables create_lex_objective(
        operations_research::sat::CpModelBuilder &cp_model,
        const int current_priority)
    {
        VectorBoolVariables objective_vars;

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
        {
            const std::vector<operations_research::sat::BoolVar> &pair_z_terms = this->all_pairs_z_terms.at(i);
            std::vector<std::vector<operations_research::sat::BoolVar>> this_phase_combos = generate_combinations(
                pair_z_terms, current_priority);

            for (const std::vector<operations_research::sat::BoolVar> &vi : this_phase_combos)
            {
                std::string var_name = "";

                for (const operations_research::sat::BoolVar &i : vi)
                {
                    var_name += i.Name();
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName("mult_" + var_name);

                // ----------------------------
                // mult ⇒ vi  (for all i)
                for (auto &v : vi)
                {
                    cp_model.AddImplication(mult, v);
                }
                // (AND(vi)) ⇒ mult  encoded as a single clause: (¬v1 ∨ ¬v2 ∨ ... ∨ mult)
                std::vector<operations_research::sat::BoolVar> clause;
                clause.reserve(vi.size() + 1);
                clause.push_back(mult);
                for (auto &v : vi)
                {
                    clause.push_back(v.Not());
                }
                cp_model.AddBoolOr(clause);
                // -----------------------------

                // If you have a fragment of length 5 (mult = 1), then there are two fragments of length 4
                // included within it. These both will be = 1.
                // This nests further (like a Russian nesting doll): if you have a fragment of length 5
                // then there are three fragments of length 3 included within. These will all be = 1, and so on.

                // cp_model.AddLessOrEqual(LinearExpr::Sum(this_phase_combos), 1).OnlyEnforceIf(this->all_obj_mults[current_priority + 1][i].Not());

                objective_vars.push_back(mult);
            }
        }

        return objective_vars;
    }

    // todo avoid recreating all priorities every lex phase - save them in base?
    void create_forward_chain_lex_objectives_as_constraints(
        operations_research::sat::CpModelBuilder &cp_model,
        const int current_priority)
    {
        int num_obj_mults = this->all_obj_mults.size();
        for (int i = current_priority; i < num_obj_mults; ++i)
        {
            std::cout << i << " " << this->all_obj_vals[i] << std::endl;

            operations_research::sat::LinearExpr freeze_expr;

            for (operations_research::sat::BoolVar &v : this->all_obj_mults[i])
            {
                // get the cloned var by proto index:
                freeze_expr += cp_model.GetBoolVarFromProtoIndex(v.index());
            }

            cp_model.AddLessOrEqual(freeze_expr, this->all_obj_vals[i]);
        }
    }

    void idk(
        operations_research::sat::CpModelBuilder &cp_model,
        const int current_priority)
    {
        int num_obj_mults = this->all_obj_mults.size();
        for (int i = num_obj_mults - 1; i > current_priority; --i)
        {
            if (this->all_obj_vals[i] == 0) continue;

            std::cout << "IDk for " << i << std::endl;

            std::vector<operations_research::sat::BoolVar> larger_fragment_obj_vars;
            for (operations_research::sat::BoolVar &v : this->all_obj_mults[i])
            {
                // get the cloned var by proto index:
                larger_fragment_obj_vars.push_back(cp_model.GetBoolVarFromProtoIndex(v.index()));
            }

            std::vector<operations_research::sat::BoolVar> smaller_fragment_obj_vars;
            for (operations_research::sat::BoolVar &v : this->all_obj_mults[i - 1])
            {
                // get the cloned var by proto index:
                smaller_fragment_obj_vars.push_back(cp_model.GetBoolVarFromProtoIndex(v.index()));
            }

            int n = this->sequence_vars_list.size();
            int num_pairs = larger_fragment_obj_vars.size() / (n * (n - 1) / 2);

            for (int j = 0; j < larger_fragment_obj_vars.size(); ++j)
            {
                int idx = j;

                if (j % (num_pairs) == 0)
                {
                    std::cout << "yes\n";
                    ++idx;
                }

                std::cout << larger_fragment_obj_vars[j].Name() << " " <<
                smaller_fragment_obj_vars[idx].Name() << " and " <<
                smaller_fragment_obj_vars[idx + 1].Name() << std::endl;

                if (j > 100) exit(1);
                // cp_model.AddImplication(larger_fragment_obj_vars[j], smaller_fragment_obj_vars[j]);
                // cp_model.AddImplication(larger_fragment_obj_vars[j], smaller_fragment_obj_vars[j + 1]);
                // cp_model.AddLessOrEqual(operations_research::sat::LinearExpr::Sum(std::vector<operations_research::sat::BoolVar>(smaller_fragment_obj_vars.begin() + j, smaller_fragment_obj_vars.begin() + j + 2)), 1).OnlyEnforceIf(larger_fragment_obj_vars[j].Not());
            }
        }
    }
};

int main(int argc, char *argv[])
{
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

        for (int i = 3; i < argc; ++i)
        {
            std::string arg = argv[i];
            std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
            if (arg == "--warmstart=true")
            {
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
    SIRIUSConfig config(
        /*log_search_progress=*/false,
        /*num_workers=*/32,
        /*linearization_level=*/3,
        /*relative_gap_limit=*/0);
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

    sirius_solver.init_cp_model(use_warmstart, sirius_solver.max_priority);

    std::cout << BLUE << "> " << RESET << "Max fragment length: " << sirius_solver.max_priority << "\n";
    std::cout << BLUE << "> " << RESET << "Current fragment length: " << current_priority << "\n";

    std::queue<char> int_vars;

    if (use_warmstart)
    {
        operations_research::sat::CpModelBuilder clone = sirius_solver.clone_cp_model(sirius_solver.base_cp_model);
        sirius_solver.set_all_var_vals_to_zero();
        sirius_solver.create_forward_chain_lex_objectives_as_constraints(clone, current_priority);
        sirius_solver.assign_var_values_from_solution();
        sirius_solver.add_hints(clone);
        sirius_solver.set_minimize_objective_value(clone, current_priority - 1);
        sirius_solver.solve_model(clone);
        sirius_solver.store_solution(current_priority - 1);

        --current_priority;
    }

    for (current_priority; current_priority > 0; --current_priority)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Setting up fragments of length " << current_priority << "...\n";

        operations_research::sat::CpModelBuilder clone = sirius_solver.clone_cp_model(sirius_solver.base_cp_model);
        // sirius_solver.create_lex_objective(clone, current_priority - 1);
        sirius_solver.create_forward_chain_lex_objectives_as_constraints(clone, current_priority);
        sirius_solver.idk(clone, current_priority);
        sirius_solver.add_hints(clone);
        sirius_solver.set_minimize_objective_value(clone, current_priority - 1);
        sirius_solver.solve_model(clone);
        sirius_solver.store_solution(current_priority - 1);
    }

    // -----

    // Variables for each base
    // Flatten sequence_vars_list
    std::vector<operations_research::sat::BoolVar> base_vars;
    for (const auto &a : sirius_solver.sequence_vars_list)
    {
        for (const auto &b : a)
        {
            for (const auto &c : b)
            {
                for (const auto &var : c)
                {
                    base_vars.push_back(var);
                }
            }
        }
    }

    // Clear
    std::queue<char>().swap(int_vars);
    for (const auto &var : base_vars)
    {
        if (operations_research::sat::SolutionIntegerValue(sirius_solver.response, var))
        {
            int_vars.push(var.Name().front());
        }
    }

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
