
#include <map>
#include <queue>
#include <memory>
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

// Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;

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
};

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

std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);

std::string translate_dna_to_protein(const std::string &dna, SIRIUSTables &tables);

std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);

std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);

void print_length_counts(const std::unordered_map<int, int> &length_counts);

void check_response(const operations_research::sat::CpSolverResponse &response);

class SIRIUSSolver {
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
    VectorBoolVariables all_mults;
    VectorBoolVariables obj_mults;
    std::vector<int> prev_obj_vals;
    VectorBoolVariables constraint_obj_mults;
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
        this->max_priority = 0;
        this->decidable_protein_length = 0;
        this->dna_with_holes = "";

        this->dna_size = this->instance.target_protein.size() * 3;
        this->target_protein_length = this->instance.target_protein.size();
        
        for (const char c : this->instance.target_protein)
        {
            max_priority += this->tables.codon_lengths_table[c];
            dna_with_holes += this->tables.invariant_codon_table[c];

            if (this->tables.skip_aa.find(c) == this->tables.skip_aa.end())
            { // if c not in skip_aa
                decidable_protein += c;
                ++decidable_protein_length;
            }
        }
    }

    void init_new_model()
    {
        all_vars.clear();

        model = std::make_unique<operations_research::sat::Model>();
        cp_model = operations_research::sat::CpModelBuilder();
        model->Add(operations_research::sat::NewSatParameters(this->config.parameters));
    }

    void build_model(const int current_priority)
    {
        this->sequence_vars_list = add_base_variables();
        this->all_mults = add_codon_constraints();
        this->all_pairs_y_terms = create_y_chained_vars();
        this->all_pairs_z_terms = create_z_chained_vars();
        this->obj_mults = create_objective(current_priority);
        this->constraint_obj_mults = create_constraint_obj_mults(current_priority);
    }

    void set_minimize_objective_value()
    {
        this->objective = operations_research::sat::LinearExpr::Sum(this->obj_mults);
        this->cp_model.Minimize(this->objective);
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

    void solve_model()
    {
        // Solve the model
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solving...\n";
        this->response = SolveCpModel(cp_model.Build(), model.get());
        check_response(this->response);

        this->prev_obj_vals.push_back(this->response.objective_value());
    }

    void store_solution()
    {
        this->map_var_name_to_val.clear();

        for (const operations_research::sat::BoolVar &v : all_vars)
        {
            this->map_var_name_to_val[v.Name()] = operations_research::sat::SolutionIntegerValue(this->response, v);
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

                    for (size_t base_idx = 0; base_idx < this->tables.codon_lengths_table.at(amino_acid); ++base_idx)
                    {
                        std::string var_name = absl::StrFormat(
                            "%c%d%d%d%d",
                            codon[base_idx],
                            sequence_n,
                            amino_acid_position,
                            codon_number,
                            base_idx);

                        operations_research::sat::BoolVar new_bool_var = this->cp_model.NewBoolVar().WithName(var_name);
                        std::cout << "Creating var " << var_name << "\n";

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
                    std::cout << "Creating var " << var_name << "\n";

                    // for (operations_research::sat::BoolVar i : codon_vars_list) {
                    //     cp_model.AddLessOrEqual(mult, i);
                    // }
                    // cp_model.AddGreaterOrEqual(mult, operations_research::sat::LinearExpr::Sum(codon_vars_list) - (codon_vars_list.size() - 1));
                    operations_research::sat::Constraint c = cp_model.AddEquality(
                        operations_research::sat::LinearExpr::Sum(codon_vars_list), codon_vars_list.size() * mult);
                    
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
                                    std::cout << "Creating var " << var_name << "\n";

                                    cp_model.AddEquality(z, y).OnlyEnforceIf(x);
                                    cp_model.AddEquality(z, 0).OnlyEnforceIf(~x);

                                    // Add z to y_terms
                                    y_terms += z;
                                    all_vars.push_back(z);
                                }
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
                operations_research::sat::BoolVar z_j = cp_model.NewBoolVar().WithName(var_name);
                std::cout << "Creating var " << var_name << "\n";

                cp_model.AddEquality(z_j, exp);

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
                std::cout << "Creating var " << var_name << "\n";

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

    VectorBoolVariables create_constraint_obj_mults(const int current_priority)
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

int main(int argc, char *argv[])
{
    // I do feel like NO ONE CAN SAVE ME I AM SO ALONE
    // And yet I CRIED, I called for HELP, forsaken
    // But now I know the only way is UNDERSTAND THE LIVING
    // OBEY THE RULE OF LIGHT and face the fear inside out
    // LOST, I found there a stone erected in line
    // WITH ONE OF THE BRIGHTEST STARS OF ALL THE NIGHT SKY VAULT
    // And I took my time, took off the moss
    // Washed away the dust and gave A NEW LEASE OF LIFE
    // Its MYSTICAL FORCE
    // I grab it now and praise this lord of earth and stone
    // Make passage for souls awaken
    // So it returns to WHERE IT'S ALWAYS BEEN with the gods

    // Systematische Identifikation Redundanter, Identisch Uebersetzter Sequenzen
    ::google::InitGoogleLogging("SIRIUS");

    // Target Protein Sequence
    // const std::string target_protein = "MM";
    // const std::string target_protein = "VTGGASK";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLIN";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLINNKDGSVLGRGHNMRFQKGSAT";
    // const std::string target_protein = "MVTGGMASKWDQKGMDIAYEEAALGYKEGGVPIGGCLINNKDGSVLGRGHNMRFQKGSATLHGEISTLENCGRLEGKVYKDTTLYTTLSPCDMCTGAIIMYGIPRCVVGENVNFKSKGEKYLQTRGHEVVVVDDERCKKIMKQFIDERPQDWFEDIGE**";
    // const std::string target_can1    = "MTNSKEDADIEEKHMYNEPVTTLFHDVEASQTHHRRGSIPLKDEKSKELYPLRSFPTRVNGEDTFSMEDGIGDEDEGEVQNAEVKRELKQRHIGMIALGGTIGTGLFIGLSTPLTNAGPVGALISYLFMGSLAYSVTQSLGEMATFIPVTSSFTVFSQRFLSPAFGAANGYMYWFSWAITFALELSVVGQVIQFWTYKVPLAAWISIFWVIITIMNLFPVKYYGEFEFWVASIKVLAIIGFLIYCFCMVCGAGVTGPVGFRYWRNPGAWGPGIISKDKNEGRFLGWVSSLINAAFTFQGTELVGITAGEAANPRKSVPRAIKKVVFRILTFYIGSLLFIGLLVPYNDPKLTQSTSYVSTSPFIIAIENSGTKVLPHIFNAVILTTIISAANSNIYVGSRILFGLSKNKLAPKFLSRTTKGGVPYIAVFVTAAFGALAYMETSTGGDKVFEWLLNITGVAGFFAWLFISISHIRFMQALKYRGISRDELPFKAKLMPGLAYYAATFMTIIIIIQGFTAFAPKFNGVSFAAAYISIFLFLAVWILFQCIFRCRFIWKIGDVDIDSDRRDIEAIVWEDHEPKTFWDKFWNVVA*";

    int num_sequences;
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

    SIRIUSTables tables;
    SIRIUSConfig config(true, 50);
    SIRIUSInstance instance(num_sequences, init_target_protein);
    SIRIUSSolver sirius_solver(instance, config, tables);
    int max_priority = sirius_solver.max_priority;

    std::queue<char> int_vars;

    for (int current_priority = max_priority; current_priority > 0; --current_priority)
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