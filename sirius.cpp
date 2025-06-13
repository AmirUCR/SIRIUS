/****************************************************************/
/*                                                              */
/*         LOST, I found there a stone erected in line          */
/*  WITH ONE OF THE BRIGHTEST STARS OF ALL THE NIGHT SKY VAULT  */
/*            And I took my time, took off the moss             */
/*      Washed away the dust and gave A NEW LEASE OF LIFE       */
/*                      Its MYSTICAL FORCE                      */
/*    I grab it now and praise this lord of earth and stone     */
/*                Make passage for souls awaken                 */
/*    So it returns to WHERE IT'S ALWAYS BEEN with the gods     */
/*                                                              */
/*                     From the Sky - Gojira                    */
/****************************************************************/

#include <map>
#include <queue>
#include <tuple>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <filesystem>

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

// TYPE DEFINITIONS
using VectorBoolVariables = std::vector<operations_research::sat::BoolVar>;
using VectorLinearExprs = std::vector<operations_research::sat::LinearExpr>;
using VectorColumnYVars = std::vector<VectorLinearExprs>;
using VectorColumnZVars = std::vector<VectorBoolVariables>;
using VectorCodonBoolVariables = std::vector<VectorBoolVariables>;
using VectorAminoacidCodonBoolVariables = std::vector<VectorCodonBoolVariables>;
using VectorSequenceAminoacidCodonBoolVariables = std::vector<VectorAminoacidCodonBoolVariables>;
// Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;
// ----------------

// CLASS FORWARD DECLARATIONS
struct SIRIUSTables;
class SIRIUSInstance;
class SIRIUSConfig;
class SIRIUSSolver;
// --------------------------

// FUNCTION DECLARATIONS
bool validate_user_prot_input(const std::string& protein, const std::unordered_map<char, std::vector<std::string>>& reduced_codon_table);
bool validate_user_num_seq_input(int num_sequences);
void print_inputs(const std::string& protein, int num_sequences);
void suppress_cout_if_quiet();
std::string elapsed_since_start();
std::string create_output_folder(const std::string& base = "sirius_out");
std::string timestamped_filename(const std::string& prefix);
std::string generate_unique_filename(const std::string& base_name);
void write_sequences_to_file_and_console(const std::vector<std::string>& sequences,
                                         const std::string& base_filename = "sirius_out.txt");
void validate_translated_proteins(const std::vector<std::string>& sequences,
                                  const std::string& target_protein,
                                  const SIRIUSTables& tables);
std::string translate_dna_to_protein(const std::string &dna, const SIRIUSTables &tables);
std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);
std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);
std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);
void print_length_counts(const std::unordered_map<int, int>& length_counts, std::ostream* file_out = nullptr);
void check_response(const operations_research::sat::CpSolverResponse &response);
// ---------------------

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
        {'G', {"A", "C", "G", "T"}},                     // {"GGT", "GGC", "GGA", "GGG"}
        {'H', {"C", "T"}},                                 // {"CAT", "CAC"}
        {'I', {"A", "C", "T"}},                            // {"ATT", "ATC", "ATA"}
        {'K', {"A", "G"}},                                 // {"AAA", "AAG"}
        {'L', {"CA", "CC", "CG", "CT", "TA", "TG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', {""}},                                       // {"ATG"}
        {'N', {"C", "T"}},                                 // {"AAT", "AAC"}
        {'P', {"A", "C", "G", "T"}},                        // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', {"A", "G"}},                                 // {"CAA", "CAG"}
        {'R', {"CT", "CC", "CA", "CG", "AA", "AG"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', {"A", "C", "G", "T"}},                     // {"ACT", "ACC", "ACA", "ACG"}
        {'V', {"A", "C", "G", "T"}},                     // {"GTT", "GTC", "GTA", "GTG"}},
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
    int decidable_protein_length;
    SIRIUSTables tables;
    std::string dna_with_holes;
    std::string target_protein;
    std::string decidable_protein;

    SIRIUSInstance(
        int n,
        std::string target_protein,
        SIRIUSTables tables)
        : n(n), target_protein(std::move(target_protein)), tables(std::move(tables))
    {
        this->dna_with_holes = "";
        this->decidable_protein = "";
        this->decidable_protein_length = 0;

        for (const char c : this->target_protein)
        {
            dna_with_holes += this->tables.invariant_codon_table[c];

            if (this->tables.skip_aa.find(c) == this->tables.skip_aa.end())
            {   // if c not in skip_aa
                decidable_protein += c;
                ++decidable_protein_length;
            }
        }
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
    int decidable_protein_length;
    std::string dna_with_holes;
    std::string decidable_protein;
    VectorBoolVariables all_vars;
    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;
    std::vector<operations_research::sat::IntVar> additive_obj_mults;
    VectorColumnYVars all_pairs_y_terms;
    VectorColumnZVars all_pairs_z_terms;
    VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;
    operations_research::sat::LinearExpr objective;

    SIRIUSSolver(SIRIUSInstance instance, SIRIUSConfig config, SIRIUSTables tables)
        : instance(std::move(instance)),
          config(std::move(config)),
          tables(std::move(tables)),
          model(std::make_unique<operations_research::sat::Model>()),
          cp_model()
    {
        this->codon_len = 3;
        this->dna_with_holes = "";
        this->decidable_protein_length = 0;

        this->dna_size = this->instance.target_protein.size() * 3;
        this->dna_with_holes = this->instance.dna_with_holes;
        this->decidable_protein = this->instance.decidable_protein;
        this->decidable_protein_length = this->instance.decidable_protein_length;
    }

    void init_new_model()
    {
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Initializing a new model...\n";
        all_vars.clear();

        // TODO model can be cloned
        // https://github.com/google/or-tools/blob/stable/ortools/sat/docs/model.md#introduction
        cp_model = operations_research::sat::CpModelBuilder();
        model = std::make_unique<operations_research::sat::Model>();
        model->Add(operations_research::sat::NewSatParameters(this->config.parameters));
    }

    void build_model()
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Construct a model...\n";
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating base vars...\n";
        this->sequence_vars_list = add_base_variables();

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating codon constraints...\n";
        this->sequence_codons_list = add_codon_constraints();

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating symmetry breaking constraints...\n";
        // add_codon_mult_anti_symmetry_constraints();
        add_codon_mult_relaxed_anti_symmetry_constraints();

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating Y chained vars...\n";
        this->all_pairs_y_terms = create_y_chained_vars();

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating Z chained vars...\n";
        this->all_pairs_z_terms = create_z_chained_vars();

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Creating additive objective...\n";
        this->additive_obj_mults = create_additive_objective();
    }

    void set_minimize_objective_value()
    {
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Setting minimize objective...\n";
        this->objective = operations_research::sat::LinearExpr::Sum(this->additive_obj_mults);
        this->cp_model.Minimize(this->objective);
    }

    void solve_model()
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solve an enigma...\n";
        this->response = SolveCpModel(cp_model.Build(), model.get());
        check_response(this->response);
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

    // std::vector<int> generate_sequence(int n_sequences, int num_codons) {
    //     std::vector<int> patterns;

    //     int n_full_cycles = n_sequences / num_codons;
    //     int n_remainder_elements = n_sequences % num_codons;
        
    //     if (n_full_cycles == 0)
    //     {
    //         std::vector<int> temp;
            
    //         for (int i = 1; i < num_codons + 1; ++i)
    //         {
    //             temp.push_back(i);
    //         }
            
    //         for (int i = num_codons - n_remainder_elements; i < num_codons; ++i)
    //         {
    //             patterns.push_back(temp.at(i));
    //         }
    //     }
    //     else
    //     {
    //         for (int i = 0; i < n_full_cycles; ++i)
    //         {
    //             for (int j = 1; j < num_codons + 1; ++j)
    //             {
    //                 patterns.push_back(j);
    //             }
    //         }
            
    //         for (int i = num_codons - n_remainder_elements; i < num_codons; ++i)
    //         {
    //             patterns.push_back(patterns.at(i));
    //         }
    //     }

    //     return patterns;
    // }

    // Function is the same as the commented above, optimized version here courtesy of ChatGPT
    std::vector<int> generate_sequence(int n_sequences, int num_codons) {
        std::vector<int> patterns;

        int n_full_cycles = n_sequences / num_codons;
        int n_remainder_elements = n_sequences % num_codons;

        // Add full cycles of 1..num_codons
        for (int i = 0; i < n_full_cycles * num_codons; ++i) {
            patterns.push_back(i % num_codons + 1);
        }

        // Handle remainder by reusing the tail from the last full cycle (or base sequence if none)
        int offset = num_codons - n_remainder_elements;
        for (int i = offset; i < num_codons; ++i) {
            patterns.push_back(i + 1);
        }

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

                std::vector<int> pattern_end_ranges = generate_sequence(n_sequences, num_codons);
    
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

                this->cp_model.AddEquality(z_j, exp);

                this->all_vars.push_back(z_j);
                this_pair_z_terms.push_back(z_j);
                colcount += 1;
            }

            all_pairs_z_terms.push_back(this_pair_z_terms);
        }

        return all_pairs_z_terms;
    }

    std::vector<operations_research::sat::IntVar> create_additive_objective()
    {
        std::vector<operations_research::sat::IntVar> factored_terms;

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i) {
            operations_research::sat::IntVar y = this->cp_model.NewIntVar({0, 1}).WithName(this->all_pairs_z_terms.at(i).front().Name());
            this->cp_model.AddEquality(y, this->all_pairs_z_terms.at(i).front());
            factored_terms.push_back(y);

            for (int j = 1; j < this->all_pairs_z_terms.at(i).size(); ++j) {
                operations_research::sat::BoolVar x = this->all_pairs_z_terms.at(i).at(j);
                operations_research::sat::IntVar z = this->cp_model.NewIntVar({0, j + 1}).WithName(x.Name() + "(1 + " + y.Name() + ")");
                
                this->cp_model.AddEquality(z, 1 + y).OnlyEnforceIf(x);
                this->cp_model.AddEquality(z, 0).OnlyEnforceIf(x.Not());
                
                y = z;
                factored_terms.push_back(y);
            }
        }

        return factored_terms;
    }
};

int main(int argc, char *argv[])
{
    // Systematische Identifikation Redundanter, Identisch Uebersetzter Sequenzen
    ::google::InitGoogleLogging("SIRIUS");
    std::cout << BLUE << "> " << RESET << "Fly to " << BLUE << "SIRIUS" << RESET << std::endl;

    bool quiet;
    int num_sequences;
    std::string init_target_protein;
    SIRIUSTables tables;

    if (argc >= 3) 
    {
        bool valid = true;
        init_target_protein = argv[1];
        try
        {
            num_sequences = std::stoi(argv[2]);
        } catch (const std::exception& e)
        {
            std::cout << RED << "> " << RESET << "Error: Number of sequences must be an integer.\n";
            valid = false;
        }

        if (!validate_user_prot_input(init_target_protein, tables.reduced_codon_table) || !validate_user_num_seq_input(num_sequences))
        {
            valid = false;
        }

        if (!valid)
        {
            std::cout << RED << "> " << RESET << "Return to Earth.\n";
            return 1;
        }

        for (int i = 3; i < argc; ++i)
        {
            std::string arg = argv[i];
            std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
            if (arg == "--quiet=true")
            {
                quiet = true;
            }
        }
    } else {
        std::cout << BLUE << "> " << RESET << "Gather your protein: ";
        std::cin >> init_target_protein;

        if (!validate_user_prot_input(init_target_protein, tables.reduced_codon_table))
        {
            std::cout << RED << "> " << RESET << "Return to Earth.\n";
            return 1;
        }

        std::cout << BLUE << "> " << RESET << "And the number of sequences: ";
        std::string num_input;
        std::cin >> num_input;

        try {
            num_sequences = std::stoi(num_input);
        } catch (const std::exception& e) {
            std::cout << RED << "> " << RESET << "Error: Number of sequences must be an integer.\n";
            std::cout << RED << "> " << RESET << "Return to Earth.\n";
            return 1;
        }

        if (!validate_user_num_seq_input(num_sequences))
        {
            std::cout << RED << "> " << RESET << "Return to Earth.\n";
            return 1;
        }
    }

    if (quiet) {
        suppress_cout_if_quiet();  // silence std::cout from here on
    }

    print_inputs(init_target_protein, num_sequences);
    
    SIRIUSConfig config(false, 16, 0, 0);
    SIRIUSInstance instance(num_sequences, init_target_protein, tables);
    SIRIUSSolver sirius_solver(instance, config, tables);

    std::queue<char> int_vars;

    // -----
    // GO!
    sirius_solver.init_new_model();
    sirius_solver.build_model();
    sirius_solver.set_minimize_objective_value();
    sirius_solver.solve_model();
    // -----

    // Variables for each base
    // Flatten sequence_vars_list
    std::vector<operations_research::sat::BoolVar> base_vars;
    for (const auto &a : sirius_solver.sequence_vars_list)
        for (const auto &b : a)
            for (const auto &c : b)
                for (const auto &var : c)
                    base_vars.push_back(var);

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

    std::string output_folder = create_output_folder();
    std::string sequences_filename = output_folder + "/" + timestamped_filename("sequences");
    std::string length_counts_filename = output_folder + "/" + timestamped_filename("length_counts");

    write_sequences_to_file_and_console(all_out_seqs, sequences_filename);
    validate_translated_proteins(all_out_seqs, init_target_protein, tables);

    auto [all_stretches, length_counts] = find_all_homologous_stretches_and_count_lengths(all_out_seqs);

    std::ofstream out_lengths(length_counts_filename);
    print_length_counts(length_counts, &out_lengths);
    out_lengths.close();

    return 0;
}

// ===============================================

bool validate_user_prot_input(const std::string& protein, const std::unordered_map<char, std::vector<std::string>>& reduced_codon_table)
{
    // Extract valid amino acid keys from the codon table
    std::unordered_set<char> valid_amino_acids;
    for (const auto& entry : reduced_codon_table)
    {
        valid_amino_acids.insert(entry.first);
    }

    // Check all characters in the protein string
    for (char c : protein)
    {
        if (!valid_amino_acids.count(c))
        {
            std::cout << RED << "> " << RESET << "Error: Invalid amino acid '" << c << "' in protein sequence.\n";
            return false;
        }
    }

    return true;
}

bool validate_user_num_seq_input(int num_sequences)
{
    // Check bounds on number of sequences
    if (num_sequences <= 1)
    {
        std::cout << RED << "> " << RESET << "Error: Must have be able to generate at least 2 sequences.\n";
        return false;
    }

    return true;
}

void print_inputs(const std::string& protein, int num_sequences)
{
    if (protein.size() > 43)  // 20 start + 3 dots + 20 end = 43 total
    {
        std::string to_print;
        to_print += protein.substr(0, 20);
        to_print += "...";
        to_print += protein.substr(protein.size() - 20);

        std::cout << BLUE << "> " << RESET << "Then in that star's heart forge " << num_sequences << "x " << to_print << std::endl;
    }
    else
    {
        std::cout << BLUE << "> " << RESET << "Then in that star's heart forge " << num_sequences << "x " << protein << std::endl;
    }
}

void suppress_cout_if_quiet() {
    static std::ofstream null_stream("/dev/null"); // use "nul" on Windows
    static std::streambuf* original_cout_buf = std::cout.rdbuf();
    std::cout.rdbuf(null_stream.rdbuf());
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

std::string create_output_folder(const std::string& base) {
    namespace fs = std::filesystem;
    std::string folder = base;
    int index = 1;
    while (fs::exists(folder)) {
        folder = base + "_" + std::to_string(index++);
    }
    fs::create_directory(folder);
    return folder;
}

std::string timestamped_filename(const std::string& prefix) {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << prefix << "_"
       << std::put_time(std::localtime(&now_c), "%Y-%m-%d_%H-%M-%S") << ".txt";
    return ss.str();
}

std::string generate_unique_filename(const std::string& base_name) {
    std::string name = base_name;
    int counter = 2;

    while (std::filesystem::exists(name)) {
        size_t dot_pos = base_name.find('.');
        if (dot_pos == std::string::npos) {
            name = base_name + "_" + std::to_string(counter);
        } else {
            name = base_name.substr(0, dot_pos) + "_" + std::to_string(counter) + base_name.substr(dot_pos);
        }
        ++counter;
    }

    return name;
}

void write_sequences_to_file_and_console(const std::vector<std::string>& sequences, const std::string& base_filename) {
    std::string filename = generate_unique_filename(base_filename);
    std::ofstream out_file(filename);
    
    if (!out_file) {
        std::cout << RED << "> " << RESET << "[" << elapsed_since_start() << "] Error: Failed to open output file.\n";
        return;
    }

    bool too_large = false;
    if (sequences.at(0).size() > 100)
    {
        too_large = true;
        std::cout << BLUE << "> " << RESET << "Resulting sequences too long to display here.\n";
    }

    for (const auto& seq : sequences) {
        if (!too_large) 
        {
            std::cout << BLUE << "> " << RESET << seq << std::endl;
        }
        
        out_file << seq << "\n";
    }

    out_file.close();
    std::cout << BLUE << "> " << RESET << "Find your sequences in " << filename << std::endl;
}

// Function to validate protein translation of sequences
void validate_translated_proteins(const std::vector<std::string>& sequences,
                                  const std::string& target_protein,
                                  const SIRIUSTables& tables) {
    for (size_t i = 0; i < sequences.size(); ++i) {
        std::string protein = translate_dna_to_protein(sequences[i], tables);

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Protein " << i + 1 << ": " << protein << "\n";
        if (protein != target_protein) {
            std::cout << RED << "> " << RESET << "Dev error: SEQ " << i << " " << sequences[i] << " NOT THE SAME\n";
        }
    }
}

std::string translate_dna_to_protein(const std::string &dna, const SIRIUSTables &tables)
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

void print_length_counts(const std::unordered_map<int, int>& length_counts, std::ostream* file_out)
{
    if (length_counts.empty()) {
        return;  // Do not print anything if there are no counts
    }
    
    std::vector<std::pair<int, int>> sorted_counts(length_counts.begin(), length_counts.end());

    std::sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto& a, const auto& b) {
                  return a.first > b.first;
              });

    std::cout << BLUE << "> " << RESET << "Fragment length counts:\n";
    for (const auto& [length, count] : sorted_counts) {
        std::cout << BLUE << "> " << RESET << "Length " << length << ": " << count << " occurrences\n";
        if (file_out) {
            *file_out << "Length " << length << ": " << count << " occurrences\n";
        }
    }
}
// -------------------------------------------------

// Print out solver response message
void check_response(const operations_research::sat::CpSolverResponse &response)
{
    if (response.status() == operations_research::sat::CpSolverStatus::OPTIMAL)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Found optimal solution.\n";
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Objective value: " << response.objective_value() << std::endl;
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::FEASIBLE)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Found feasible solution.\n";
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Objective value: " << response.objective_value() << std::endl;
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE)
    {
        std::cout << RED << "> DEV ERROR " << RESET << "[" << elapsed_since_start() << "] Infeasible.\n";
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID)
    {
        std::cout << RED << "> DEV ERROR " << RESET << "[" << elapsed_since_start() << "] MODEL INVALID.\n";
    }
}
