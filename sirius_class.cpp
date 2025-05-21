
#include <map>
#include <queue>
#include <deque>
#include <tuple>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <unordered_map>

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

struct SIRIUSTables {
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
        {'V', "_T_"}, // {"GTT", "GTC", "GTA", "GTG"}
        {'W', "TGG"}, // TGG
        {'Y', "TA_"}, // {"TAT", "TAC"}
        {'*', "T__"}  // {"TAA", "TAG", "TGA"}
    };

    std::set<char> skip_aa = {'M', 'W'};

    std::unordered_map<char, std::vector<std::string>> codon_table = {
        {'A', {"T", "C", "A", "G"}},                       // {"GCT", "GCC", "GCA", "GCG"}
        {'C', {"T", "C"}},                                 // {"TGT", "TGC"}
        {'D', {"T", "C"}},                                 // {"GAT", "GAC"}
        {'E', {"A", "G"}},                                 // {"GAA", "GAG"}
        {'F', {"T", "C"}},                                 // {"TTT", "TTC"}
        {'G', {"T", "C", "A", "G"}},                       // {"GGT", "GGC", "GGA", "GGG"}
        {'H', {"T", "C"}},                                 // {"CAT", "CAC"}
        {'I', {"T", "C", "A"}},                            // {"ATT", "ATC", "ATA"}
        {'K', {"A", "G"}},                                 // {"AAA", "AAG"}
        {'L', {"TA", "TG", "CT", "CC", "CA", "CG"}},       // {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
        {'M', {""}},                                       // {"ATG"}
        {'N', {"T", "C"}},                                 // {"AAT", "AAC"}
        {'P', {"T", "C", "A", "G"}},                       // {"CCT", "CCC", "CCA", "CCG"}
        {'Q', {"A", "G"}},                                 // {"CAA", "CAG"}
        {'R', {"CT", "CC", "CA", "CG", "AA", "AG"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', {"T", "C", "A", "G"}},                       // {"ACT", "ACC", "ACA", "ACG"}
        {'V', {"GT", "GC", "GA", "GG"}},                   // {"GTT", "GTC", "GTA", "GTG"}},
        {'W', {""}},                                       // {"TGG"}
        {'Y', {"T", "C"}},                                 // {"TAT", "TAC"}},
        {'*', {"AA", "AG", "GA"}}                          // {"TAA", "TAG", "TGA"}
    };

    std::unordered_map<char, int> codon_lengths_table = {
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
        {'V', 2},
        {'W', 0},
        {'Y', 1},
        {'*', 2}};

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
}

class SIRIUSInstance {
    public:
    int n;
    std::string target_protein;

    SIRIUSInstance(int n, std::string target_protein)
        : n(n), target_protein(std::move(target_protein))
    {
        // Pass
    }
};

class SIRIUSConfig {
    public:
    operations_research::sat::SatParameters parameters;

    SIRIUSConfig(bool log_search_progress, unsigned int num_workers)
    {
        parameters.set_log_search_progress(log_search_progress);
        parameters.set_num_workers(num_workers);
    }
};

class SIRIUSSolver {
    public:
    SIRIUSTables tables;
    SIRIUSConfig config;
    SIRIUSInstance instance;
    operations_research::sat::Model model;
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
    VectorBoolVariables all_mults;
    VectorBoolVariables obj_mults;
    std::vector<int> prev_obj_vals;
    VectorBoolVariables constraint_obj_mults;
    VectorColumnYVars all_pairs_y_terms;
    VectorColumnZVars all_pairs_z_terms;
    VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;
    operations_research::sat::LinearExpr objective;

    SIRIUSSolver(SIRIUSInstance instance, SIRIUSConfig config, SIRIUSTables tables) {
        this->codon_len = 3;
        this->max_priority = 0;
        this->decideable_protein_length = 0;
        this->dna_with_holes = "";

        this->tables = tables;
        this->config = config;
        this->instance = instance;
        this->dna_size = instance.target_protein.size() * 3;
        this->target_protein_length = instance.target_protein.size();
        
        for (const char c : instance.target_protein)
        {
            max_priority += tables.codon_lengths_table[c];
            dna_with_holes += tables.invariant_codon_table[c];

            if (tables.skip_aa.find(c) == tables.skip_aa.end())
            { // if c not in skip_aa
                decidable_protein += c;
                ++decidable_protein_length;
            }
        }

        // this->build_model();
    }

    private:
    void build_model(const int current_priority)
    {
        this->sequence_vars_list = add_base_variables();
        this->all_mults = add_codon_constraints();
        this->all_pairs_y_terms = create_y_chained_vars();
        this->all_pairs_z_terms = create_z_chained_vars();
        this->obj_mults = create_objective(current_priority);
        this->constraint_obj_mults = constraint_obj_mults(current_priority);
    }

    void set_minimize_objective_value()
    {
        this->objective = operations_research::sat::LinearExpr::Sum(this->obj_mults);
        this->cp_model.Minimize(this->objective);
    }

    void solve_model()
    {

    }

    void add_hints()
    {
        for (const operations_research::sat::BoolVar &v : all_vars)
        {
            if (this->map_var_name_to_val.find(v.Name()) != this->map_var_name_to_val.end())
            {
                this->cp_model.AddHint(v, this->map_var_name_to_val.at(v.Name()));
            }
        }
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

                for (size_t codon_number = 0; codon_number < this->tables.codon_table.at(amino_acid).size(); ++codon_number)
                {
                    const std::string &codon = this->tables.codon_table.at(amino_acid).at(codon_number);
                   VectorBoolVariables base_vars_list;

                    for (size_t base_idx = 0; base_idx < codon_lengths_table.at(amino_acid); ++base_idx)
                    {
                        std::string var_name = absl::StrFormat(
                            "%c%d%d%d%d",
                            codon[base_idx],
                            sequence_n,
                            amino_acid_position,
                            codon_number,
                            base_idx);

                        operations_research::sat::BoolVar new_bool_var = cp_model.NewBoolVar().WithName(var_name);

                        all_vars.push_back(new_bool_var);
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

    VectorBoolVariables add_codon_constraints()
    {
        VectorBoolVariables all_mults;

        // Constraints for Valid Codons
        for (const auto &this_sequence_vars_list_of_list_it : this->sequence_vars_list)
        {
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

                    operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                    // std::cout << "Creating var " << var_name << "\n";

                    // for (operations_research::sat::BoolVar i : codon_vars_list) {
                    //     cp_model.AddLessOrEqual(mult, i);
                    // }
                    // cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(codon_vars_list) - (codon_vars_list.size() - 1));
                    operations_research::sat::Constraint c = cp_model.AddEquality(
                        operations_research::sat::LinearExpr::Sum(codon_vars_list), codon_vars_list.size() * mult);

                    // if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                    // {
                    //     cp_model.AddHint(mult, map_var_name_to_val.at(var_name));
                    //     std::cout << "Adding hint for mult " << var_name << "\n";
                    // }
                    
                    all_vars.push_back(mult);
                    codon_mult_list.push_back(mult);
                    all_mults.push_back(mult);
                }

                if (codon_mult_list.size() != 0)
                {
                    // Ensure only one codon auxiliary variable is chosen
                    operations_research::sat::Constraint c = cp_model.AddEquality(
                        operations_research::sat::LinearExpr::Sum(codon_mult_list), 1);
                }
            }
        }

        return all_mults;
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
                    for (int codon_i = 0; codon_i < this->tables.codon_lengths_table.at(this->decidable_protein.at(aa_pos_i)); ++codon_i)
                    {
                        operations_research::sat::LinearExpr y_terms;

                        for (size_t codon_position_i = 0; codon_position_i < this->tables.codon_table.at(this->decidable_protein.at(aa_pos_i)).size(); ++codon_position_i)
                        {
                            operations_research::sat::BoolVar x = this->sequence_vars_list.at(s).at(aa_pos_i).at(codon_position_i).at(codon_i);

                            for (size_t codon_position_j = 0; codon_position_j < this->tables.codon_table.at(this->decidable_protein.at(aa_pos_i)).size(); ++codon_position_j)
                            {
                                operations_research::sat::BoolVar y = this->sequence_vars_list.at(t).at(aa_pos_i).at(codon_position_j).at(codon_i);

                                if (x.Name().front() == y.Name().front())
                                {
                                    std::string var_name = x.Name() + y.Name();

                                    operations_research::sat::BoolVar z = cp_model.NewBoolVar().WithName(var_name);
                                    // std::cout << "Creating var " << var_name << "\n";

                                    cp_model.AddEquality(z, y).OnlyEnforceIf(x);
                                    cp_model.AddEquality(z, 0).OnlyEnforceIf(~x);

                                    // if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                                    // {
                                    //     cp_model.AddHint(z, map_var_name_to_val.at(var_name));
                                    //     std::cout << "Adding hint for z " << var_name << "\n";
                                    // }

                                    // Add z to y_terms
                                    y_terms += z;
                                    all_vars.push_back(z);
                                }
                            }
                        }
                        this_pairs_y_terms.push_back(y_terms);
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
                operations_research::sat::BoolVar z_j = cp_model.NewBoolVar().WithName(var_name);
                cp_model.AddEquality(z_j, y_terms);

                all_vars.push_back(z_j);
                this_pair_z_terms.push_back(z_j);
                colcount += 1;
            }

            all_pairs_z_terms.push_back(this_pair_z_terms);
        }

        return all_pairs_z_terms;
    }

    VectorBoolVariables create_objective(const int current_priority)
    {
        VectorBoolVariables all_mults;

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
        {
            std::vector<std::vector<operations_research::sat::BoolVar>> combos = generate_combinations(
                this->all_pairs_z_terms.at(i), current_priority);

            for (const std::vector<operations_research::sat::BoolVar> vi : combos)
            {
                std::string var_name = "";

                for (const operations_research::sat::BoolVar i : vi)
                {
                    var_name += i.Name();
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);

                for (const operations_research::sat::BoolVar i : vi)
                {
                    cp_model.AddLessOrEqual(mult, i);
                }
                cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));

                all_vars.push_back(mult);
                all_mults.push_back(mult);
            }
        }

        return all_mults;
    }

    VectorBoolVariables constraint_obj_mults(const int current_priority)
    {
        int counter = 0;
        VectorBoolVariables constraint_obj_mults;

        for (int go_down = this->max_priority; go_down > current_priority; --go_down)
        {
            VectorBoolVariables obj_mults;

            for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
            {
                std::vector<VectorBoolVariables> combos = generate_combinations(
                    this->all_pairs_z_terms.at(i), go_down);

                for (VectorBoolVariables vi : combos)
                {
                    std::string var_name = "";

                    for (const operations_research::sat::BoolVar &i : vi)
                    {
                        var_name += i.Name();
                    }

                    operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                    std::cout << "Creating var " << var_name << "\n";

                    for (const operations_research::sat::BoolVar &i : vi)
                    {
                        cp_model.AddLessOrEqual(mult, i);
                    }
                    cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));

                    // if (this->map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                    // {
                    //     cp_model.AddHint(mult, map_var_name_to_val.at(var_name));
                    //     std::cout << "Adding hint for obj_mult " << var_name << "\n";
                    // }

                    all_vars.push_back(mult);
                    obj_mults.push_back(mult);
                    constraint_obj_mults.push_back(mult);
                }
            }

            if (this->prev_obj_vals.size() > counter)
            {
                cp_model.AddLessOrEqual(operations_research::sat::LinearExpr::Sum(obj_mults), this->prev_obj_vals.at(counter));
                counter++;
            }
        }

        return constraint_obj_mults;
    }
};



std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);

std::string translate_dna_to_protein(const std::string &dna);

// Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;

std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);

std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);

void print_length_counts(const std::unordered_map<int, int> &length_counts);

void check_response(const operations_research::sat::CpSolverResponse &response);

std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> create_variables(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::string target_protein,
    const int prot_len,
    const int num_sequences,
    const std::unordered_map<char, int> &codon_lengths_table,
    const std::unordered_map<char, std::vector<std::string>> &codon_table);

std::vector<operations_research::sat::BoolVar> create_codon_constraints(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> &sequence_vars_list);

std::vector<std::vector<operations_research::sat::BoolVar>> create_z_chained_vars(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    std::vector<operations_research::sat::BoolVar> &z_vars,
    const std::string target_protein,
    const int prot_len,
    const int num_sequences,
    const std::unordered_map<char, int> &codon_lengths_table,
    const std::unordered_map<char, std::vector<std::string>> &codon_table,
    const std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> &sequence_vars_list);

std::vector<operations_research::sat::BoolVar> create_objective(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<operations_research::sat::BoolVar>> &all_pairs_z_terms,
    const int current_priority);

std::vector<operations_research::sat::BoolVar> create_constraints_from_objectives(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<operations_research::sat::BoolVar>> &all_pairs_z_terms,
    const std::vector<int> &prev_obj_vals,
    const int current_priority,
    const int max_priority);

int main(int argc, char *argv[])
{
    // Now from a solar system to another I fly!
    // Systematische Identifikation Redundanter, Identisch Uebersetzter Sequenzen
    ::google::InitGoogleLogging("TO SIRIUS");

    // Target Protein Sequence
    // const std::string target_protein = "MM";
    // const std::string target_protein = "VTGGASK";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLIN";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLINNKDGSVLGRGHNMRFQKGSAT";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLINNKDGSVLGRGHNMRFQKGSATLHGEISTLENCGRLEGKVYKDTTLYTTLSPCDMCTGAIIMYGIPRCVVGENVNFKSKGEKYLQTRGHEVVVVDDERCKKIMKQFIDERPQDWFEDIGE**";
    // const std::string target_can1    = "MTNSKEDADIEEKHMYNEPVTTLFHDVEASQTHHRRGSIPLKDEKSKELYPLRSFPTRVNGEDTFSMEDGIGDEDEGEVQNAEVKRELKQRHIGMIALGGTIGTGLFIGLSTPLTNAGPVGALISYLFMGSLAYSVTQSLGEMATFIPVTSSFTVFSQRFLSPAFGAANGYMYWFSWAITFALELSVVGQVIQFWTYKVPLAAWISIFWVIITIMNLFPVKYYGEFEFWVASIKVLAIIGFLIYCFCMVCGAGVTGPVGFRYWRNPGAWGPGIISKDKNEGRFLGWVSSLINAAFTFQGTELVGITAGEAANPRKSVPRAIKKVVFRILTFYIGSLLFIGLLVPYNDPKLTQSTSYVSTSPFIIAIENSGTKVLPHIFNAVILTTIISAANSNIYVGSRILFGLSKNKLAPKFLSRTTKGGVPYIAVFVTAAFGALAYMETSTGGDKVFEWLLNITGVAGFFAWLFISISHIRFMQALKYRGISRDELPFKAKLMPGLAYYAATFMTIIIIIQGFTAFAPKFNGVSFAAAYISIFLFLAVWILFQCIFRCRFIWKIGDVDIDSDRRDIEAIVWEDHEPKTFWDKFWNVVA*";

    int num_sequences;
    int max_priority = 0;
    std::string init_target_protein;

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
        // Use command-line arguments
        init_target_protein = argv[1];
        num_sequences = std::stoi(argv[2]);
    }

    std::cout << BLUE << "> " << RESET << "Creating " << num_sequences << "x " << init_target_protein << std::endl;

    const int CODON_LEN = 3;
    const int INIT_PROT_LEN = init_target_protein.size();

    std::string dna_with_holes = "";
    for (auto i : init_target_protein)
    {
        dna_with_holes += invariant_codon_table[i];
        max_priority += codon_lengths_table[i];
    }

    int PROT_LEN = 0;
    std::string target_protein;
    for (char c : init_target_protein)
    {
        if (skip_aa.find(c) == skip_aa.end())
        { // if c not in skip_aa
            target_protein += c;
            ++PROT_LEN;
        }
    }

    std::queue<char> int_vars;
    std::vector<int> prev_obj_vals;
    std::unordered_map<std::string, int> map_var_name_to_val;

    operations_research::sat::SatParameters parameters;
    parameters.set_log_search_progress(true);
    parameters.set_num_workers(50);

    for (int current_priority = max_priority; current_priority > 0; --current_priority)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Setting up fragments of length " << current_priority << "...\n";

        // Create a model
        operations_research::sat::Model model;
        operations_research::sat::CpModelBuilder cp_model;
        operations_research::sat::CpSolverResponse response;

        model.Add(operations_research::sat::NewSatParameters(parameters));

        // Variables for each base
        std::vector<
            std::vector<
                std::vector<
                    std::vector<
                        operations_research::sat::BoolVar>>>>
            sequence_vars_list = create_variables(cp_model,
                                                  map_var_name_to_val,
                                                  target_protein,
                                                  PROT_LEN,
                                                  num_sequences,
                                                  codon_lengths_table,
                                                  codon_table);

        std::vector<operations_research::sat::BoolVar> all_mults = create_codon_constraints(cp_model,
                                                                                            map_var_name_to_val,
                                                                                            sequence_vars_list);

        std::vector<operations_research::sat::BoolVar> z_vars;
        std::vector<
            std::vector<
                operations_research::sat::BoolVar>>
            all_pairs_z_terms = create_z_chained_vars(cp_model,
                                                      map_var_name_to_val,
                                                      z_vars,
                                                      target_protein,
                                                      PROT_LEN,
                                                      num_sequences,
                                                      codon_lengths_table,
                                                      codon_table,
                                                      sequence_vars_list);

        std::vector<operations_research::sat::BoolVar> obj_mults = create_objective(cp_model,
                                                                                    map_var_name_to_val,
                                                                                    all_pairs_z_terms,
                                                                                    current_priority);

        operations_research::sat::LinearExpr objective = operations_research::sat::LinearExpr::Sum(obj_mults);

        std::vector<operations_research::sat::BoolVar> constraint_obj_mults = create_constraints_from_objectives(cp_model,
                                                                                                                 map_var_name_to_val,
                                                                                                                 all_pairs_z_terms,
                                                                                                                 prev_obj_vals,
                                                                                                                 current_priority,
                                                                                                                 max_priority);

        cp_model.Minimize(objective);

        // Solve the model
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solving...\n";
        response = SolveCpModel(cp_model.Build(), &model);
        check_response(response);

        prev_obj_vals.push_back(response.objective_value());

        // Flatten sequence_vars_list
        std::vector<operations_research::sat::BoolVar> all_vars;
        for (const auto &a : sequence_vars_list)
            for (const auto &b : a)
                for (const auto &c : b)
                    for (const auto &var : c)
                        all_vars.push_back(var);

        // Store var values to be hinted on next solve
        map_var_name_to_val.clear();
        for (const operations_research::sat::BoolVar &var : all_vars)
        {
            map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
        }

        for (const operations_research::sat::BoolVar &var : all_mults)
        {
            map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
        }

        for (const std::vector<operations_research::sat::BoolVar> &vec : all_pairs_z_terms)
        {
            for (const operations_research::sat::BoolVar &var : vec)
            {
                map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
            }
        }

        for (const operations_research::sat::BoolVar &var : obj_mults)
        {
            map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
        }

        for (const operations_research::sat::BoolVar &var : constraint_obj_mults)
        {
            map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
        }

        for (const operations_research::sat::BoolVar &var : z_vars)
        {
            map_var_name_to_val[var.Name()] = operations_research::sat::SolutionIntegerValue(response, var);
        }

        // Clear
        std::queue<char>().swap(int_vars);
        for (const auto &var : all_vars)
        {
            if (operations_research::sat::SolutionIntegerValue(response, var))
            {
                int_vars.push(var.Name().front());
                // std::cout << var.Name() << std::endl;
            }
        }
    }

    // -----

    int counter = 0;
    std::string seq;

    std::vector<std::string> all_out_seqs;
    for (int seq_n = 0; seq_n < num_sequences; ++seq_n)
    {
        std::string this_seq = "";
        for (int i = 0; i < INIT_PROT_LEN * CODON_LEN; ++i)
        {
            if (dna_with_holes[i] == '_')
            {
                this_seq += int_vars.front();
                int_vars.pop();
            }
            else
            {
                this_seq += dna_with_holes[i];
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
        std::string protein = translate_dna_to_protein(all_out_seqs[i]);

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

std::string translate_dna_to_protein(const std::string &dna)
{
    std::string protein;
    for (size_t i = 0; i + 2 < dna.size(); i += 3)
    {
        std::string codon = dna.substr(i, 3);
        auto it = translate_codon_table.find(codon);
        if (it != translate_codon_table.end())
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

// Create variables for each base and assign to cp_model
std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> create_variables(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::string target_protein,
    const int prot_len,
    const int num_sequences,
    const std::unordered_map<char, int> &codon_lengths_table,
    const std::unordered_map<char, std::vector<std::string>> &codon_table)
{

    // Variables for each base
    std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> sequence_vars_list;

    for (int sequence_n = 0; sequence_n < num_sequences; ++sequence_n)
    {
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> this_sequence_vars_list_of_list;

        for (int amino_acid_position = 0; amino_acid_position < prot_len; ++amino_acid_position)
        {
            char amino_acid = target_protein[amino_acid_position];
            std::vector<std::vector<operations_research::sat::BoolVar>> codon_vars_list;

            for (size_t codon_number = 0; codon_number < codon_table.at(amino_acid).size(); ++codon_number)
            {
                const std::string &codon = codon_table.at(amino_acid).at(codon_number);
                std::vector<operations_research::sat::BoolVar> base_vars_list;

                for (size_t base_idx = 0; base_idx < codon_lengths_table.at(amino_acid); ++base_idx)
                {
                    std::string var_name = absl::StrFormat(
                        "%c%d%d%d%d",
                        codon[base_idx],
                        sequence_n,
                        amino_acid_position,
                        codon_number,
                        base_idx);

                    operations_research::sat::BoolVar new_bool_var = cp_model.NewBoolVar().WithName(var_name);
                    std::cout << "Creating var " << var_name << "\n";

                    if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                    {
                        cp_model.AddHint(new_bool_var, map_var_name_to_val.at(var_name));
                        std::cout << "Adding hint for base var " << var_name << "\n";
                    }

                    base_vars_list.push_back(new_bool_var);
                }
                if (base_vars_list.size() > 0)
                {
                    codon_vars_list.push_back(base_vars_list);
                }
            }
            if (codon_vars_list.size() > 0)
            {
                this_sequence_vars_list_of_list.push_back(codon_vars_list);
            }
        }
        if (this_sequence_vars_list_of_list.size() > 0)
        {
            sequence_vars_list.push_back(this_sequence_vars_list_of_list);
        }
    }

    return sequence_vars_list;
}

// Create auxiliary variables such that an aux variable represents up to 3 base variables
// which the solver has to consider. If fewer bases than 3 in a codon must be considered, chain
// those together instead of all 3.
// Also constrain the model such that only one codon (or fewer considered bases) must be chosen
// per amino acid.
std::vector<operations_research::sat::BoolVar> create_codon_constraints(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> &sequence_vars_list)
{
    std::vector<operations_research::sat::BoolVar> all_mults;

    // Constraints for Valid Codons
    for (const auto &this_sequence_vars_list_of_list_it : sequence_vars_list)
    {
        for (const auto &codon_vars_list_of_lists : this_sequence_vars_list_of_list_it)
        {
            std::vector<operations_research::sat::BoolVar> codon_mult_list;

            for (const auto &codon_vars_list : codon_vars_list_of_lists)
            {
                std::string var_name;
                for (const auto &var : codon_vars_list)
                {
                    var_name += var.Name();
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                std::cout << "Creating var " << var_name << "\n";

                // for (operations_research::sat::BoolVar i : codon_vars_list) {
                //     cp_model.AddLessOrEqual(mult, i);
                // }
                // cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(codon_vars_list) - (codon_vars_list.size() - 1));
                operations_research::sat::Constraint c = cp_model.AddEquality(
                    operations_research::sat::LinearExpr::Sum(codon_vars_list), codon_vars_list.size() * mult);

                if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                {
                    cp_model.AddHint(mult, map_var_name_to_val.at(var_name));
                    std::cout << "Adding hint for mult " << var_name << "\n";
                }

                codon_mult_list.push_back(mult);
                all_mults.push_back(mult);
            }

            if (codon_mult_list.size() != 0)
            {
                // Ensure only one codon auxiliary variable is chosen
                operations_research::sat::Constraint c = cp_model.AddEquality(
                    operations_research::sat::LinearExpr::Sum(codon_mult_list), 1);
            }
        }
    }

    return all_mults;
}

// Create z variables. A z variable is a sum of all base-multiplied aux variables in a column of bases.
// For example,
// S1 = CAA
//      CAG
// S2 = CAA
//      CAG
// Then z1 = (CC + CC + CC + CC) -- This will not exist actually because of invariant_codon_table
//                                     and the fact that we mask columns of bases with no choice to make.
//      z2 = (AA + AA + AA + AA) -- This will not exist in the solver.
//      z3 = (AA + GG) -- this will exist because a decision has to be made.
std::vector<std::vector<operations_research::sat::BoolVar>> create_z_chained_vars(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    std::vector<operations_research::sat::BoolVar> &z_vars,
    const std::string target_protein,
    const int prot_len,
    const int num_sequences,
    const std::unordered_map<char, int> &codon_lengths_table,
    const std::unordered_map<char, std::vector<std::string>> &codon_table,
    const std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> &sequence_vars_list)
{
    std::vector<std::vector<operations_research::sat::BoolVar>> all_pairs_z_terms;

    for (int s = 0; s < num_sequences; ++s)
    {
        for (int t = s + 1; t < num_sequences; ++t)
        {
            std::vector<operations_research::sat::BoolVar> this_pair_z_terms;
            int colcount = 1;

            for (int aa_pos_i = 0; aa_pos_i < prot_len; ++aa_pos_i)
            {
                for (int codon_i = 0; codon_i < codon_lengths_table.at(target_protein.at(aa_pos_i)); ++codon_i)
                {
                    operations_research::sat::LinearExpr y_terms;

                    for (size_t codon_position_i = 0; codon_position_i < codon_table.at(target_protein.at(aa_pos_i)).size(); ++codon_position_i)
                    {
                        operations_research::sat::BoolVar x = sequence_vars_list.at(s).at(aa_pos_i).at(codon_position_i).at(codon_i);

                        for (size_t codon_position_j = 0; codon_position_j < codon_table.at(target_protein.at(aa_pos_i)).size(); ++codon_position_j)
                        {
                            operations_research::sat::BoolVar y = sequence_vars_list.at(t).at(aa_pos_i).at(codon_position_j).at(codon_i);

                            if (x.Name().front() == y.Name().front())
                            {
                                std::string var_name = x.Name() + y.Name();

                                operations_research::sat::BoolVar z = cp_model.NewBoolVar().WithName(var_name);
                                std::cout << "Creating var " << var_name << "\n";

                                cp_model.AddEquality(z, y).OnlyEnforceIf(x);
                                cp_model.AddEquality(z, 0).OnlyEnforceIf(~x);

                                if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                                {
                                    cp_model.AddHint(z, map_var_name_to_val.at(var_name));
                                    std::cout << "Adding hint for z " << var_name << "\n";
                                }

                                // Add z to y_terms
                                y_terms += z;
                                z_vars.push_back(z);
                            }
                        }
                    }

                    std::string var_name = absl::StrFormat("z%d%d%d", s, t, colcount);
                    operations_research::sat::BoolVar z_j = cp_model.NewBoolVar().WithName(var_name);
                    std::cout << "Creating var " << var_name << "\n";
                    cp_model.AddEquality(z_j, y_terms);

                    if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                    {
                        cp_model.AddHint(z_j, map_var_name_to_val.at(var_name));
                        std::cout << "Adding hint for z_j " << var_name << "\n";
                    }

                    this_pair_z_terms.push_back(z_j);
                    colcount += 1;
                }
            }

            all_pairs_z_terms.push_back(this_pair_z_terms);
        }
    }

    return all_pairs_z_terms;
}

// Chain together z terms into auxiliary variables depending on the length of the objective priority.
// For example, when current_priority == 10, then mult == z1 * z2 * z3 * z4 * ... * z10
std::vector<operations_research::sat::BoolVar> create_objective(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<operations_research::sat::BoolVar>> &all_pairs_z_terms,
    const int current_priority)
{
    std::vector<operations_research::sat::BoolVar> all_mults;

    for (int i = 0; i < all_pairs_z_terms.size(); ++i)
    {
        std::vector<std::vector<operations_research::sat::BoolVar>> combos = generate_combinations(
            all_pairs_z_terms.at(i), current_priority);

        for (std::vector<operations_research::sat::BoolVar> vi : combos)
        {
            std::string var_name = "";

            for (auto i : vi)
            {
                var_name += i.Name();
            }

            operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
            std::cout << "Creating var " << var_name << "\n";

            for (operations_research::sat::BoolVar i : vi)
            {
                cp_model.AddLessOrEqual(mult, i);
            }
            cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));

            if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
            {
                cp_model.AddHint(mult, map_var_name_to_val.at(var_name));
                std::cout << "Adding hint for obj_mult " << var_name << "\n";
            }

            all_mults.push_back(mult);
        }
    }

    return all_mults;
}

// Chain together z terms into auxiliary variables depending on the length of the objective priority,
// but in a forward fashion such that previously optimized objectives show up here again as
// constraints with their old values in prev_obj_vals.
std::vector<operations_research::sat::BoolVar> create_constraints_from_objectives(
    operations_research::sat::CpModelBuilder &cp_model,
    const std::unordered_map<std::string, int> &map_var_name_to_val,
    const std::vector<std::vector<operations_research::sat::BoolVar>> &all_pairs_z_terms,
    const std::vector<int> &prev_obj_vals,
    const int current_priority,
    const int max_priority)
{
    int counter = 0;
    int slack = 0;
    std::vector<operations_research::sat::BoolVar> constraint_obj_mults;

    for (int go_down = max_priority; go_down > current_priority; --go_down)
    {
        std::vector<operations_research::sat::BoolVar> obj_mults;

        for (int i = 0; i < all_pairs_z_terms.size(); ++i)
        {
            std::vector<std::vector<operations_research::sat::BoolVar>> combos = generate_combinations(
                all_pairs_z_terms.at(i), go_down);

            for (std::vector<operations_research::sat::BoolVar> vi : combos)
            {
                std::string var_name = "";

                for (const operations_research::sat::BoolVar &i : vi)
                {
                    var_name += i.Name();
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName(var_name);
                std::cout << "Creating var " << var_name << "\n";

                for (const operations_research::sat::BoolVar &i : vi)
                {
                    cp_model.AddLessOrEqual(mult, i);
                }
                cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(vi) - (vi.size() - 1));

                if (map_var_name_to_val.find(var_name) != map_var_name_to_val.end())
                {
                    cp_model.AddHint(mult, map_var_name_to_val.at(var_name));
                    std::cout << "Adding hint for obj_mult " << var_name << "\n";
                }

                obj_mults.push_back(mult);
                constraint_obj_mults.push_back(mult);
            }
        }

        cp_model.AddLessOrEqual(operations_research::sat::LinearExpr::Sum(obj_mults), prev_obj_vals.at(counter));
        counter++;
    }

    return constraint_obj_mults;
}
