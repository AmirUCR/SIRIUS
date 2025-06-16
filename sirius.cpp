/****************************************************************/
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
#include <cmath>
#include <queue>
#include <tuple>
#include <cctype>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <utility>
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

bool quiet = false;

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
class SIRIUSTables;
class SIRIUSInstance;
class SIRIUSConfig;
class SIRIUSSolver;
// --------------------------

// FUNCTION DECLARATIONS
void scan_for_quiet_flag(int argc, char *argv[]);
void print_info(const std::string &message);
void print_info_newline(const std::string &message);
void print_warning(const std::string &message);
void print_warning_newline(const std::string &message);
void print_error(const std::string &message);
void print_error_newline(const std::string &message);
std::string format_error(const std::string &message);
[[noreturn]] void throw_formatted_error(const std::string &message);

void validate_user_prot_input(const std::string &protein,
                              const std::unordered_map<char, std::vector<std::string>> &reduced_codon_table);
void validate_user_num_seq_input(int num_sequences);
void print_inputs(const std::string &protein, int num_sequences);

// SIRIUSInstance gather_inputs_from_yaml(const std::string& path, const SIRIUSTables& tables);
std::pair<SIRIUSInstance, SIRIUSConfig> gather_inputs_from_flags(int argc, char *argv[], const SIRIUSTables &tables);
std::pair<SIRIUSInstance, SIRIUSConfig> gather_inputs_interactively(SIRIUSTables &tables);

std::string elapsed_since_start();
std::string create_output_folder(const std::string &base = "sirius_out");
std::string timestamped_filename(const std::string &prefix);
std::string generate_unique_filename(const std::string &base_name);
void write_sequences_to_file_and_console(const std::vector<std::string> &sequences,
                                         const std::string &base_filename = "sirius_out.txt");
void validate_translated_proteins(const std::vector<std::string> &sequences,
                                  const std::string &target_protein,
                                  const SIRIUSTables &tables);
std::string translate_dna_to_protein(const std::string &dna, const SIRIUSTables &tables);

std::vector<std::vector<operations_research::sat::BoolVar>> generate_combinations(
    const std::vector<operations_research::sat::BoolVar> &vars,
    const int combination_size);
std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);

std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);
void print_length_counts(const std::unordered_map<int, int> &length_counts, std::ostream *file_out = nullptr);
void check_response(const operations_research::sat::CpSolverResponse &response);
// ---------------------

class SIRIUSTables
{
public:
    std::unordered_map<char, std::unordered_map<std::string, double>> rscu_map;

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

    std::set<char> skip_aa = {'M', 'W'};

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
        {'R', {"CT", "CC", "CA", "CG", "AA", "AG"}},       // {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
        {'S', {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}}, // {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
        {'T', {"A", "C", "G", "T"}},                       // {"ACT", "ACC", "ACA", "ACG"}
        {'V', {"A", "C", "G", "T"}},                       // {"GTT", "GTC", "GTA", "GTG"}},
        {'W', {""}},                                       // {"TGG"}
        {'Y', {"C", "T"}},                                 // {"TAT", "TAC"}},
        {'*', {"AA", "AG", "GA"}}                          // {"TAA", "TAG", "TGA"}
    };

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

    void build_rscu_map_from_csv(const std::string &filename)
    {
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            throw_formatted_error("Error: Cannot open RSCU CSV file.");
        }

        std::string line;
        std::getline(infile, line); // Skip header
        int line_number = 1;

        while (std::getline(infile, line))
        {
            ++line_number;

            std::stringstream ss(line);
            std::string aa_str, tmp1, codon, tmp2, rscu_str, tmp3;

            std::getline(ss, aa_str, ',');   // AmOneLet
            std::getline(ss, tmp1, ',');     // AmAcid
            std::getline(ss, codon, ',');    // Codon
            std::getline(ss, tmp2, ',');     // Percentage
            std::getline(ss, rscu_str, ','); // RSCU
            std::getline(ss, tmp3, ',');     // GC3

            if (aa_str.empty() || codon.size() != 3 || rscu_str.empty())
            {
                throw_formatted_error("Error: Malformed RSCU CSV line " + std::to_string(line_number) + ": " + line);
            }

            char aa = aa_str[0];
            double rscu = 0.0;
            try
            {
                rscu = std::stod(rscu_str);
            }
            catch (const std::invalid_argument &)
            {
                throw_formatted_error("Error: RSCU value on line " + std::to_string(line_number) + " could not be converted.");
            }

            auto it = this->invariant_codon_table.find(aa);
            if (it == this->invariant_codon_table.end())
            {
                throw_formatted_error("Error: Amino acid '" + std::to_string(aa) + "' on line " + std::to_string(line_number) + " is unknown.");
            }

            const std::string &pattern = it->second;
            std::string varying_part;

            for (size_t i = 0; i < pattern.size(); ++i)
            {
                if (pattern[i] == '_')
                {
                    varying_part += codon[i];
                }
            }

            if (!varying_part.empty())
            {
                this->rscu_map[aa][varying_part] = rscu;
            }
        }
    }
};

class SIRIUSInstance
{
public:
    int n;
    int decidable_protein_length;
    double hard_rscu_threshold;
    double soft_rscu_threshold;
    double gc_ending_rscu_threshold;
    double rscu_alpha;
    double max_low_rscu_ratio;
    bool hard_filter_by_rscu = false;
    bool soft_filter_by_rscu = false;
    ;

    std::string codon_usage_path;

    SIRIUSTables tables;
    std::string dna_with_holes;
    std::string decidable_protein;
    std::string init_target_protein;

    SIRIUSInstance() = default;

    SIRIUSInstance(
        int n,
        std::string init_target_protein,
        SIRIUSTables tables,
        std::string codon_usage_path,
        double hard_rscu_threshold,
        double soft_rscu_threshold,
        double gc_ending_rscu_threshold,
        double rscu_alpha,
        double max_low_rscu_ratio)
        : n(n),
          init_target_protein(std::move(init_target_protein)),
          tables(std::move(tables)),
          codon_usage_path(std::move(codon_usage_path)),
          hard_rscu_threshold(hard_rscu_threshold),
          soft_rscu_threshold(soft_rscu_threshold),
          gc_ending_rscu_threshold(gc_ending_rscu_threshold),
          rscu_alpha(rscu_alpha),
          max_low_rscu_ratio(max_low_rscu_ratio)
    {
        std::transform(
            this->init_target_protein.begin(),
            this->init_target_protein.end(),
            this->init_target_protein.begin(), ::toupper);

        if (this->hard_rscu_threshold > 0)
        {
            this->hard_filter_by_rscu = true;
        }
        if (this->soft_rscu_threshold > 0)
        {
            this->soft_filter_by_rscu = true;
        }

        this->dna_with_holes = "";
        this->decidable_protein = "";
        this->decidable_protein_length = 0;

        for (const char c : this->init_target_protein)
        {
            dna_with_holes += this->tables.invariant_codon_table[c];

            if (this->tables.skip_aa.find(c) == this->tables.skip_aa.end())
            { // if c not in skip_aa
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
    SIRIUSConfig() = default;

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
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

    SIRIUSTables tables;
    SIRIUSConfig config;
    SIRIUSInstance instance;

    operations_research::sat::CpModelBuilder cp_model;
    operations_research::sat::CpSolverResponse response;
    std::unique_ptr<operations_research::sat::Model> model;

    int dna_size;
    int codon_len;
    int decidable_protein_length;

    std::string dna_with_holes;
    std::string decidable_protein;

    VectorBoolVariables all_vars;
    VectorColumnYVars all_pairs_y_terms;
    VectorColumnZVars all_pairs_z_terms;
    operations_research::sat::LinearExpr objective;
    VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;
    std::vector<operations_research::sat::IntVar> additive_obj_mults;
    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;

    SIRIUSSolver(SIRIUSInstance instance, SIRIUSConfig config, SIRIUSTables tables)
        : instance(std::move(instance)),
          config(std::move(config)),
          tables(std::move(tables)),
          model(std::make_unique<operations_research::sat::Model>()),
          cp_model(),
          gen(std::random_device{}()), // gen(42), // todo make injectable
          dis(0.0, 1.0)
    {
        this->codon_len = 3;
        this->dna_with_holes = "";
        this->decidable_protein_length = 0;

        this->dna_size = this->instance.init_target_protein.size() * 3;
        this->dna_with_holes = this->instance.dna_with_holes;
        this->decidable_protein = this->instance.decidable_protein;
        this->decidable_protein_length = this->instance.decidable_protein_length;
    }

    void init_new_model()
    {
        cp_model = operations_research::sat::CpModelBuilder();
        model = std::make_unique<operations_research::sat::Model>();
        model->Add(operations_research::sat::NewSatParameters(this->config.parameters));
    }

    void build_model()
    {
        print_info_newline("[" + elapsed_since_start() + "] Construct a model...");
        // print_info_newline("[" + elapsed_since_start() + "] Creating base vars...");
        this->sequence_vars_list = (this->instance.hard_filter_by_rscu || this->instance.soft_filter_by_rscu)
                                       ? add_base_variables_from_rscu()
                                       : add_base_variables_from_prot();

        // print_info_newline("[" + elapsed_since_start() + "] Creating codon constraints...");
        this->sequence_codons_list = add_codon_constraints();

        // print_info_newline("[" + elapsed_since_start() + "] Creating symmetry breaking constraints...");
        if (!this->instance.soft_filter_by_rscu)
        {
            add_codon_mult_relaxed_anti_symmetry_constraints();
        }

        // print_info_newline("[" + elapsed_since_start() + "] Creating Y chained vars...");
        this->all_pairs_y_terms = create_y_chained_vars();

        // print_info_newline("[" + elapsed_since_start() + "] Creating Z chained vars...");
        this->all_pairs_z_terms = create_z_chained_vars();

        // print_info_newline("[" + elapsed_since_start() + "] Creating additive objective...");
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

    std::vector<std::vector<std::string>> generate_sequence_codons_with_soft_rscu()
    {
        // Map each amino acid to all positions it occurs at in the decidable protein
        std::unordered_map<char, std::vector<int>> aa_to_positions;
        std::unordered_map<char, std::unordered_map<std::string, int>> num_issued_codons;
        std::unordered_map<char, std::unordered_map<std::string, int>> max_allowed_per_codon;
        std::vector<std::vector<std::string>> sequence_codons(this->decidable_protein_length);
        std::unordered_set<char> unique_aminoacids(this->decidable_protein.begin(), this->decidable_protein.end());

        for (const char &aa : unique_aminoacids)
        {
            int count_of_aa = std::count(this->decidable_protein.begin(), this->decidable_protein.end(), aa);

            for (const std::string &codon : this->tables.reduced_codon_table.at(aa))
            {
                double rscu = this->tables.rscu_map.at(aa).at(codon);

                num_issued_codons[aa][codon] = 0;

                // Marked for Russian roulette!
                if (rscu < this->instance.soft_rscu_threshold)
                {
                    // Allow low-RSCU codon at most X% of times it's possible for this amino acid
                    // e.g., 30% of the total positions where this AA occurs
                    int max_count = std::ceil(this->instance.max_low_rscu_ratio * count_of_aa);
                    max_allowed_per_codon[aa][codon] = max_count;
                }
            }
        }

        // Save the occurrence positions of the protein's AA's
        for (int i = 0; i < this->decidable_protein.size(); ++i)
        {
            aa_to_positions[this->decidable_protein.at(i)].push_back(i);
        }

        // Go through each aa and its positions
        for (auto &[aa, positions] : aa_to_positions)
        {
            // Random shuffle so we don't run out of
            // low-RSCU codons by the end of the protein
            std::shuffle(positions.begin(), positions.end(), this->gen);

            for (int p : positions)
            {
                std::vector<std::string> codons;

                double largest_rscu = 0;
                bool at_least_one_passed = false;
                std::string codon_w_largest_rscu = "";

                for (const std::string &codon : this->tables.reduced_codon_table.at(aa))
                {
                    double rscu = this->tables.rscu_map.at(aa).at(codon);

                    if (rscu > largest_rscu)
                    {
                        largest_rscu = rscu;
                        codon_w_largest_rscu = codon;
                    }

                    // Russian roulette time!
                    if (rscu < this->instance.soft_rscu_threshold)
                    {
                        // Enforce max use cap
                        if (num_issued_codons.at(aa).at(codon) >= max_allowed_per_codon.at(aa).at(codon))
                        {
                            continue;
                        }

                        double normalized = std::max(0.0, std::min(1.0, rscu / this->instance.soft_rscu_threshold));
                        double probability = std::exp(-this->instance.rscu_alpha * (1.0 - normalized));

                        // Jetzt hast du Pech gehabt.
                        if (this->dis(this->gen) > probability)
                        {
                            continue;
                        }

                        num_issued_codons.at(aa).at(codon) += 1;
                    }

                    codons.push_back(codon);
                    at_least_one_passed = true;
                }

                // If nothing passes the threshold, add the codon w/ highest RSCU
                if (!at_least_one_passed)
                {
                    // std::ostringstream oss;
                    // oss << "Warning: None of the codons for amino acid "
                    //     << aa << " position " << std::to_string(p + 1) << " pass the set RSCU threshold ("
                    //     << this->instance.rscu_threshold << "). Letting "
                    //     << codon_w_largest_rscu << " with the largest RSCU ("
                    //     << largest_rscu << ") through.";
                    // print_warning_newline(oss.str());

                    codons.push_back(codon_w_largest_rscu);
                }

                sequence_codons.at(p) = codons;
            }
        }

        return sequence_codons;
    }

    std::vector<std::vector<std::string>> generate_sequence_codons_with_hard_rscu()
    {
        std::vector<std::vector<std::string>> sequence_codons(this->decidable_protein_length);

        for (int i = 0; i < this->decidable_protein.size(); ++i)
        {
            char aa = this->decidable_protein.at(i);

            std::vector<std::string> codons;

            double largest_rscu = 0;
            bool at_least_one_passed = false;
            std::string codon_w_largest_rscu = "";

            for (const std::string &codon : this->tables.reduced_codon_table.at(aa))
            {
                std::string full_codon = "";

                int reduced_codon_idx = 0;
                for (const char &c : this->tables.invariant_codon_table.at(aa))
                {
                    if (c != '_')
                    {
                        full_codon += c;
                    }
                    else
                    {
                        full_codon += codon.at(reduced_codon_idx);
                        ++reduced_codon_idx;
                    }
                }

                double rscu = this->tables.rscu_map.at(aa).at(codon);

                if (rscu > largest_rscu)
                {
                    largest_rscu = rscu;
                    codon_w_largest_rscu = codon;
                }

                // Axed.
                if (rscu < this->instance.hard_rscu_threshold)
                {
                    continue;
                }

                if (rscu < this->instance.gc_ending_rscu_threshold &&
                    (full_codon.back() != 'G' && full_codon.back() != 'C'))
                {
                    continue;
                }

                codons.push_back(codon);
                at_least_one_passed = true;
            }

            // If nothing passes the threshold, add the codon w/ highest RSCU
            if (!at_least_one_passed)
            {
                codons.push_back(codon_w_largest_rscu);
            }

            sequence_codons.push_back(codons);
        }

        return sequence_codons;
    }

    std::vector<std::vector<std::string>> generate_sequence_codons_with_rscu()
    {
        return (this->instance.soft_filter_by_rscu) ? generate_sequence_codons_with_soft_rscu() : generate_sequence_codons_with_hard_rscu();
    }

    VectorSequenceAminoacidCodonBoolVariables add_base_variables_from_rscu()
    {
        // Variables for each base
        VectorSequenceAminoacidCodonBoolVariables sequence_vars_list;

        for (int sequence_n = 0; sequence_n < this->instance.n; ++sequence_n)
        {
            VectorAminoacidCodonBoolVariables this_sequence_vars_list_of_list;
            std::vector<std::vector<std::string>> sequence_with_codons = generate_sequence_codons_with_rscu();

            for (int amino_acid_position = 0; amino_acid_position < sequence_with_codons.size(); ++amino_acid_position)
            {
                VectorCodonBoolVariables codon_vars_list;

                for (size_t codon_number = 0; codon_number < sequence_with_codons.at(amino_acid_position).size(); ++codon_number)
                {
                    const std::string &codon = sequence_with_codons.at(amino_acid_position).at(codon_number);

                    VectorBoolVariables base_vars_list;
                    for (size_t base_idx = 0; base_idx < codon.size(); ++base_idx)
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

    VectorSequenceAminoacidCodonBoolVariables add_base_variables_from_prot()
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

    std::vector<int> generate_sequence(int n_sequences, int num_codons)
    {
        std::vector<int> patterns;

        int n_full_cycles = n_sequences / num_codons;
        int n_remainder_elements = n_sequences % num_codons;

        // Add full cycles of 1..num_codons
        for (int i = 0; i < n_full_cycles * num_codons; ++i)
        {
            patterns.push_back(i % num_codons + 1);
        }

        // Handle remainder by reusing the tail from the last full cycle (or base sequence if none)
        int offset = num_codons - n_remainder_elements;
        for (int i = offset; i < num_codons; ++i)
        {
            patterns.push_back(i + 1);
        }

        return patterns;
    }

    void add_codon_mult_relaxed_anti_symmetry_constraints()
    {
        int n_sequences = this->sequence_codons_list.size();

        for (int k = 0; k < n_sequences - 1; ++k)
        {
            int seq_length = this->sequence_codons_list[k].size();
            std::vector<operations_research::sat::BoolVar> prefix_equal(seq_length);

            for (int pos = 0; pos < seq_length; ++pos)
            {
                // Determine if seq[k] and seq[k+1] codons at this position are equal
                prefix_equal[pos] = this->cp_model.NewBoolVar();
                int num_codons = this->sequence_codons_list[k][pos].size();

                std::vector<int> pattern_end_ranges = generate_sequence(n_sequences, num_codons);

                // Equality at this position across all codon variables
                std::vector<operations_research::sat::BoolVar> codon_equals;
                for (int c = 0; c < num_codons; ++c)
                {
                    operations_research::sat::BoolVar codon_match = cp_model.NewBoolVar();
                    this->cp_model.AddEquality(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k + 1][pos][c]).OnlyEnforceIf(codon_match);
                    this->cp_model.AddNotEqual(this->sequence_codons_list[k][pos][c], this->sequence_codons_list[k + 1][pos][c]).OnlyEnforceIf(codon_match.Not());
                    codon_equals.push_back(codon_match);
                }
                // All codon bits must match to declare equality at this position
                this->cp_model.AddBoolAnd(codon_equals).OnlyEnforceIf(prefix_equal[pos]);

                std::vector<operations_research::sat::BoolVar> codon_differs;
                for (const auto &eq : codon_equals)
                {
                    codon_differs.push_back(eq.Not());
                }
                this->cp_model.AddBoolOr(codon_differs).OnlyEnforceIf(prefix_equal[pos].Not());

                operations_research::sat::BoolVar lex_less_or_equal = this->cp_model.NewBoolVar();

                // Ensures prefix_equal[pos] = AND(codon_equals)
                // Enforce lex ordering at this position
                // Create ordering constraint: seq[k] <= seq[k+1]
                this->cp_model.AddImplication(prefix_equal[pos].Not(), lex_less_or_equal);
                std::vector<operations_research::sat::BoolVar> selector_vars; // one for each (i,j) pair

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

                        for (int bi = 0; bi < num_codons; ++bi)
                        {
                            if (bi != i && bi != j)
                            {
                                this->cp_model.AddLessOrEqual(this->sequence_codons_list[k][pos][bi], this->sequence_codons_list[k + 1][pos][bi])
                                    .OnlyEnforceIf({lex_less_or_equal, selector});
                            }
                        }
                    }
                }

                // Create a dummy variable to accumulate the sum of selector_vars
                operations_research::sat::LinearExpr selector_sum;
                for (const auto &var : selector_vars)
                {
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

                        for (size_t codon_position_i = 0; codon_position_i < this->sequence_vars_list.at(s).at(aa_pos_i).size(); ++codon_position_i)
                        {
                            operations_research::sat::BoolVar x = this->sequence_vars_list.at(s).at(aa_pos_i).at(codon_position_i).at(codon_i);

                            for (size_t codon_position_j = 0; codon_position_j < this->sequence_vars_list.at(t).at(aa_pos_i).size(); ++codon_position_j)
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

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
        {
            operations_research::sat::IntVar y = this->cp_model.NewIntVar({0, 1}).WithName(this->all_pairs_z_terms.at(i).front().Name());
            this->cp_model.AddEquality(y, this->all_pairs_z_terms.at(i).front());
            factored_terms.push_back(y);

            for (int j = 1; j < this->all_pairs_z_terms.at(i).size(); ++j)
            {
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
    scan_for_quiet_flag(argc, argv);
    print_info_newline(std::string("Fly to ") + BLUE + std::string("SIRIUS") + RESET);

    SIRIUSConfig config;
    SIRIUSTables tables;
    SIRIUSInstance instance;

    if (argc == 1)
    {
        // instance = gather_inputs_from_yaml("config.yaml", tables);
        std::tie(instance, config) = gather_inputs_interactively(tables); // todo
    }
    else if (argc == 2 && std::string(argv[1]) == "-i")
    {
        std::tie(instance, config) = gather_inputs_interactively(tables);
    }
    else
    {
        std::tie(instance, config) = gather_inputs_from_flags(argc, argv, tables);
    }

    if (instance.hard_filter_by_rscu || instance.soft_filter_by_rscu)
    {
        tables.build_rscu_map_from_csv(instance.codon_usage_path);
    }

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

    std::string seq;

    std::vector<std::string> all_out_seqs;
    for (int seq_n = 0; seq_n < instance.n; ++seq_n)
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
    validate_translated_proteins(all_out_seqs, instance.init_target_protein, tables);

    auto [all_stretches, length_counts] = find_all_homologous_stretches_and_count_lengths(all_out_seqs);

    try
    {
        std::ofstream out_lengths(length_counts_filename);
        print_length_counts(length_counts, &out_lengths);
        out_lengths.close();
    }
    catch (const std::exception &e)
    {
        throw e;
    }
}

// ===============================================
void scan_for_quiet_flag(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        std::transform(arg.begin(), arg.end(), arg.begin(), ::tolower);

        if (arg.find("--quiet=") == 0)
        {
            // quiet is global
            quiet = arg.substr(8) == "true";
        }
    }
}

void print_info(const std::string &message)
{
    if (!quiet)
    {
        std::cout << BLUE << "> " << RESET << message;
    }
}

void print_info_newline(const std::string &message)
{
    if (!quiet)
    {
        std::cout << BLUE << "> " << RESET << message << "\n";
    }
}

void print_warning(const std::string &message)
{
    if (!quiet)
    {
        std::cout << ORANGE << "> " << RESET << message;
    }
}

void print_warning_newline(const std::string &message)
{
    if (!quiet)
    {
        std::cout << ORANGE << "> " << RESET << message << "\n";
    }
}

void print_error(const std::string &message)
{
    if (!quiet)
    {
        std::cout << RED << "> " << RESET << message;
    }
}

void print_error_newline(const std::string &message)
{
    if (!quiet)
    {
        std::cout << RED << "> " << RESET << message << "\n";
    }
}

std::string format_error(const std::string &message)
{
    return RED + std::string("> ") + RESET + message;
}

[[noreturn]] void throw_formatted_error(const std::string &message)
{
    if (!quiet)
    {
        std::cerr << format_error(message) << "\n";
        std::cout << format_error("Back to Earth.") << "\n";
    }
    exit(1);
}

void validate_user_prot_input(
    const std::string &protein,
    const std::unordered_map<char, std::vector<std::string>> &reduced_codon_table)
{
    std::unordered_set<char> invalid_aa;

    std::string upper_prot = protein;
    std::transform(upper_prot.begin(), upper_prot.end(), upper_prot.begin(), ::toupper);

    for (char aa : upper_prot)
    {
        if (reduced_codon_table.find(aa) == reduced_codon_table.end())
        {
            invalid_aa.insert(aa);
        }
    }

    if (!invalid_aa.empty())
    {
        std::ostringstream oss;
        oss << "Invalid amino acids in protein sequence: ";
        for (char aa : invalid_aa)
        {
            oss << "'" << aa << "' ";
        }
        throw_formatted_error(oss.str());
    }
}

void validate_user_num_seq_input(int num_sequences)
{
    // Check bounds on number of sequences
    if (num_sequences <= 1)
    {
        throw_formatted_error("Must generate at least 2 sequences.");
    }
}

void validate_file_exists(const std::string &filename, const std::string &message)
{
    if (!std::filesystem::exists(filename))
    {
        throw_formatted_error("Error: file for " + message + " (\"" + filename + "\") does not exist.");
    }
}

void print_inputs(const std::string &protein, int num_sequences)
{
    // 20 start + 3 dots + 20 end = 43 total
    if (protein.size() > 43)
    {
        std::string to_print;
        to_print += protein.substr(0, 20);
        to_print += "...";
        to_print += protein.substr(protein.size() - 20);

        print_info_newline("Now, in that star's heart, forge " + std::to_string(num_sequences) + "x " + to_print);
    }
    else
    {
        print_info_newline("Now, in that star's heart, forge " + std::to_string(num_sequences) + "x " + protein);
    }
}

// SIRIUSInstance gather_inputs_from_yaml(const std::string& path, const SIRIUSTables& tables)
// {
//     YAML::Node config = YAML::LoadFile(path);

//     std::string codon_usage_path = ""; // todo

//     std::string prot = config["Prot"].as<std::string>();
//     int n = config["num_sequences"].as<int>();
//     double rscu_thresh = config["rscu_threshold"].as<double>();
//     double alpha = config["rscu_alpha"].as<double>();
//     double max_ratio = config["max_low_rscu_ratio"].as<double>();
//     bool hard = config["hard_filter_by_rscu"].as<bool>();
//     bool soft = config["soft_filter_by_rscu"].as<bool>();

//     return SIRIUSInstance(n, prot, tables, codon_usage_path, rscu_thresh, alpha, max_ratio, hard, soft);
// }

std::pair<SIRIUSInstance, SIRIUSConfig> gather_inputs_from_flags(int argc, char *argv[], const SIRIUSTables &tables)
{
    std::string just_input;

    int num_sequences = 2;
    std::string init_target_protein = "MALEEINENSTERN";

    int num_workers = 16;
    bool show_ortools_log = false;
    double relative_gap_limit = 0.0;

    bool user_set_alpha = false;
    bool user_set_rscu_ratio = false;
    bool user_set_gc_end_rscu_threshold = false;
    bool user_set_hard_rscu_threshold = false, user_set_soft_rscu_threshold = false;

    std::string codon_usage_path = "";

    double rscu_alpha = 0, max_low_rscu_ratio = 0.0;
    double hard_rscu_threshold = 0.0, soft_rscu_threshold = 0.0, gc_end_rscu_threshold = 0.0;

    if (argc == 3)
    {
        std::string first_arg = argv[1];
        std::string second_arg = argv[2];

        if (first_arg.find("--prot=") == 0)
        {
            init_target_protein = first_arg.substr(7);
        }
        else
        {
            init_target_protein = first_arg;
        }

        if (second_arg.find("--n=") == 0)
        {
            try
            {
                num_sequences = std::stoi(second_arg.substr(4));
            }
            catch (...)
            {
                throw_formatted_error("Error: Number of sequences must be a positive integer.");
            }
        }
        else
        {
            try
            {
                num_sequences = std::stoi(second_arg);
            }
            catch (...)
            {
                throw_formatted_error("Error: Number of sequences must be a positive integer.");
            }
        }

        validate_user_prot_input(init_target_protein, tables.reduced_codon_table);
        validate_user_num_seq_input(num_sequences);

        print_inputs(init_target_protein, num_sequences);

        SIRIUSConfig config(show_ortools_log, num_workers, 0, relative_gap_limit);
        SIRIUSInstance instance(num_sequences, init_target_protein, tables, "", 0, 0, 0, false, false);
        return {instance, config};
    }

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg.find("--prot=") == 0)
        {
            init_target_protein = arg.substr(7);
        }
        else if (arg.find("--n=") == 0)
        {
            try
            {
                num_sequences = std::stoi(arg.substr(4));
            }
            catch (...)
            {
                throw_formatted_error("Error: Number of sequences must be a positive integer.");
            }
        }
        else if (arg.find("--relative_gap_limit=") == 0)
        {
            try
            {
                relative_gap_limit = std::stod(arg.substr(21));
            }
            catch (...)
            {
                throw_formatted_error("Error: Relative gap limit must be a positive number.");
            }
        }
        else if (arg.find("--num_workers=") == 0)
        {
            try
            {
                num_workers = std::stoi(arg.substr(14));
            }
            catch (...)
            {
                throw_formatted_error("Error: Number of CPU workers must be a positive integer.");
            }
        }
        else if (arg.find("--show_ortools_log=") == 0)
        {
            show_ortools_log = arg.substr(19) == "true";
        }
        else if (arg.find("--soft_rscu_thresh=") == 0)
        {
            try
            {
                soft_rscu_threshold = std::stod(arg.substr(19));
                user_set_soft_rscu_threshold = true;
            }
            catch (...)
            {
                throw_formatted_error("Error: Soft RSCU threshold must be a positive number.");
            }
        }
        else if (arg.find("--hard_rscu_thresh=") == 0)
        {
            try
            {
                hard_rscu_threshold = std::stod(arg.substr(19));
                user_set_hard_rscu_threshold = true;
            }
            catch (...)
            {
                throw_formatted_error("Error: Hard RSCU threshold must be a positive number.");
            }
        }
        else if (arg.find("--gc_end_rscu_thresh=") == 0)
        {
            try
            {
                gc_end_rscu_threshold = std::stod(arg.substr(21));
                user_set_gc_end_rscu_threshold = true;
            }
            catch (...)
            {
                throw_formatted_error("Error: RSCU threshold must be a positive number.");
            }
        }
        else if (arg.find("--rscu_alpha=") == 0)
        {
            try
            {
                rscu_alpha = std::stod(arg.substr(13));
                user_set_alpha = true;
            }
            catch (...)
            {
                throw_formatted_error("Error: RSCU's alpha must be a number.");
            }
        }
        else if (arg.find("--max_low_rscu_ratio=") == 0)
        {
            try
            {
                max_low_rscu_ratio = std::stod(arg.substr(21));
                if (max_low_rscu_ratio < 0 || max_low_rscu_ratio > 1)
                {
                    throw_formatted_error("Error: Max RSCU ratio must be a number between 0-1.");
                }
                user_set_rscu_ratio = true;
            }
            catch (...)
            {
                throw_formatted_error("Error: Max RSCU ratio must be a number between 0-1.");
            }
        }
        else if (arg.find("--codon_usage_fpath=") == 0)
        {
            codon_usage_path = arg.substr(20);
            validate_file_exists(codon_usage_path, "codon usage file");
        }
        else
        {
            throw_formatted_error("Error: Unrecognized flag: " + arg);
        }
    }

    validate_user_prot_input(init_target_protein, tables.reduced_codon_table);
    validate_user_num_seq_input(num_sequences);

    if (user_set_hard_rscu_threshold && user_set_soft_rscu_threshold)
    {
        throw_formatted_error("Error: Cannot both soft- and hard-filter RSCU. Define a value for only one.");
    }

    if ((user_set_hard_rscu_threshold || user_set_soft_rscu_threshold) && codon_usage_path.empty())
    {
        throw_formatted_error("Error: You must provide a valid path for the codon usage file (--codon_usage_fpath=) when using RSCU filtering.");
    }

    if (user_set_soft_rscu_threshold && !user_set_rscu_ratio)
    {
        std::ostringstream oss;
        oss << "Warning: You did not specify an RSCU ratio for soft filtering. Using default ("
            << std::fixed << std::setprecision(1) << max_low_rscu_ratio << ").";

        print_warning_newline(oss.str());
    }

    if (user_set_soft_rscu_threshold && !user_set_alpha)
    {
        std::ostringstream oss;
        oss << "Warning: You did not specify an RSCU alpha for soft filtering. Using default ("
            << std::fixed << std::setprecision(1) << std::showpoint << rscu_alpha << ").";

        print_warning_newline(oss.str());
    }

    print_inputs(init_target_protein, num_sequences);

    SIRIUSInstance instance(
        num_sequences,
        init_target_protein,
        tables,
        codon_usage_path,
        hard_rscu_threshold,
        soft_rscu_threshold,
        gc_end_rscu_threshold,
        rscu_alpha,
        max_low_rscu_ratio);

    SIRIUSConfig config(show_ortools_log, num_workers, 0, relative_gap_limit);

    return {instance, config};
}

std::pair<SIRIUSInstance, SIRIUSConfig> gather_inputs_interactively(SIRIUSTables &tables)
{
    int num_sequences;
    std::string init_target_protein;
    std::string codon_usage_path = "";

    double rscu_alpha = 10.0;
    double max_low_rscu_ratio = 0.3;
    double hard_rscu_threshold = 0.5;
    double soft_rscu_threshold = 0.7;
    double gc_ending_rscu_threshold = 0.5;

    bool use_hard_rscu = false;
    bool use_soft_rscu = false;

    std::string input;

    print_info("Gather your protein: ");
    std::getline(std::cin, init_target_protein);
    if (init_target_protein.empty())
    {
        init_target_protein = "MALEEINENSTERN";
        print_info_newline("Setting protein to " + init_target_protein);
    }
    else
    {
        validate_user_prot_input(init_target_protein, tables.reduced_codon_table);
    }

    print_info("And the number of sequences: ");
    std::getline(std::cin, input);
    if (input.empty())
    {
        num_sequences = 2;
        print_info_newline("Setting number of sequences to " + std::to_string(num_sequences));
    }
    else
    {
        try
        {
            num_sequences = std::stoi(input);
        }
        catch (...)
        {
            throw_formatted_error("Error: Number of sequences must be an integer.");
        }
        validate_user_num_seq_input(num_sequences);
    }

    // ---- Hard filter by RSCU? ----
    print_info("Hard filter by RSCU? (yes or no) [Default no]: ");
    std::getline(std::cin, input);
    if (input.empty() || input == "no")
    {
        use_hard_rscu = false;
        hard_rscu_threshold = 0;
    }
    else if (input == "yes")
    {
        use_hard_rscu = true;
    }
    else
    {
        print_info_newline("Hä?! Bold of you to assume \"" + input + "\" for hard RSCU filter is valid. Taking that as a no.");
        use_hard_rscu = false;
        hard_rscu_threshold = 0;
    }

    if (!use_hard_rscu)
    {
        // ---- Soft filter by RSCU? ----
        print_info("Soft filter by RSCU? (yes or no) [Default no]: ");
        std::getline(std::cin, input);
        if (input.empty() || input == "no")
        {
            use_soft_rscu = false;
            soft_rscu_threshold = 0;
        }
        else if (input == "yes")
        {
            use_soft_rscu = true;
        }
        else
        {
            print_info_newline("Hä?! Bold of you to assume \"" + input + "\" for soft RSCU filter is valid. Taking that as a no.");
            use_soft_rscu = false;
            soft_rscu_threshold = 0;
        }

        if (use_soft_rscu)
        {
            print_info("Soft RSCU threshold [Default 0.7]: ");
            std::getline(std::cin, input);
            if (!input.empty())
            {
                try
                {
                    soft_rscu_threshold = std::stod(input);
                    hard_rscu_threshold = 0;
                }
                catch (...)
                {
                    throw_formatted_error("Error: Soft RSCU threshold must be a number.");
                }
            }

            print_info("Soft RSCU alpha [Default 10.0]: ");
            std::getline(std::cin, input);
            if (!input.empty())
            {
                try
                {
                    rscu_alpha = std::stod(input);
                }
                catch (...)
                {
                    throw_formatted_error("Error: Alpha must be a number.");
                }
            }

            print_info("Max low-RSCU ratio [Default 0.3]: ");
            std::getline(std::cin, input);
            if (!input.empty())
            {
                try
                {
                    max_low_rscu_ratio = std::stod(input);
                }
                catch (...)
                {
                    throw_formatted_error("Error: Max low-RSCU ratio must be a number.");
                }
            }
        }
    }
    else
    {
        print_info("Hard RSCU threshold [Default 0.5]: ");
        std::getline(std::cin, input);
        if (!input.empty())
        {
            try
            {
                hard_rscu_threshold = std::stod(input);
                soft_rscu_threshold = 0;
            }
            catch (...)
            {
                throw_formatted_error("Error: Hard RSCU threshold must be a number.");
            }
        }

        print_info("GC-ending RSCU threshold [Default 0.5]: ");
        std::getline(std::cin, input);
        if (!input.empty())
        {
            try
            {
                gc_ending_rscu_threshold = std::stod(input);
            }
            catch (...)
            {
                throw_formatted_error("Error: GC-ending RSCU threshold must be a number.");
            }
        }

        print_info("Codon usage file path: ");
        std::getline(std::cin, codon_usage_path);
        validate_file_exists(codon_usage_path, "codon usage file");
    }

    // Codon usage required if soft filter is enabled too
    if (use_soft_rscu && codon_usage_path.empty())
    {
        print_info("Codon usage file path: ");
        std::getline(std::cin, codon_usage_path);
        validate_file_exists(codon_usage_path, "codon usage file");
    }

    print_inputs(init_target_protein, num_sequences);

    SIRIUSInstance instance(
        num_sequences,
        init_target_protein,
        tables,
        codon_usage_path,
        hard_rscu_threshold,
        soft_rscu_threshold,
        gc_ending_rscu_threshold,
        rscu_alpha,
        max_low_rscu_ratio);

    SIRIUSConfig config(false, 16, 0, 0.0); // fixed defaults

    return {instance, config};
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

std::string create_output_folder(const std::string &base)
{
    namespace fs = std::filesystem;
    std::string folder = base;
    int index = 1;
    while (fs::exists(folder))
    {
        folder = base + "_" + std::to_string(index++);
    }
    fs::create_directory(folder);
    return folder;
}

std::string timestamped_filename(const std::string &prefix)
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << prefix << "_"
       << std::put_time(std::localtime(&now_c), "%Y-%m-%d_%H-%M-%S") << ".txt";
    return ss.str();
}

std::string generate_unique_filename(const std::string &base_name)
{
    std::string name = base_name;
    int counter = 2;

    while (std::filesystem::exists(name))
    {
        size_t dot_pos = base_name.find('.');
        if (dot_pos == std::string::npos)
        {
            name = base_name + "_" + std::to_string(counter);
        }
        else
        {
            name = base_name.substr(0, dot_pos) + "_" + std::to_string(counter) + base_name.substr(dot_pos);
        }
        ++counter;
    }

    return name;
}

void write_sequences_to_file_and_console(const std::vector<std::string> &sequences, const std::string &base_filename)
{
    std::string filename = generate_unique_filename(base_filename);
    std::ofstream out_file(filename);

    if (!out_file)
    {
        std::cout << RED << "> " << RESET << "[" << elapsed_since_start() << "] Error: Failed to open output file.\n";
        return;
    }

    bool too_large = false;
    if (sequences.at(0).size() > 100)
    {
        too_large = true;
        std::cout << BLUE << "> " << RESET << "Resulting sequences too long to display here.\n";
    }

    for (const auto &seq : sequences)
    {
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
void validate_translated_proteins(const std::vector<std::string> &sequences,
                                  const std::string &target_protein,
                                  const SIRIUSTables &tables)
{
    for (size_t i = 0; i < sequences.size(); ++i)
    {
        std::string protein = translate_dna_to_protein(sequences[i], tables);

        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Protein " << i + 1 << ": " << protein << "\n";
        if (protein != target_protein)
        {
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

void print_length_counts(const std::unordered_map<int, int> &length_counts, std::ostream *file_out)
{
    if (length_counts.empty())
    {
        return; // Do not print anything if there are no counts
    }

    std::vector<std::pair<int, int>> sorted_counts(length_counts.begin(), length_counts.end());

    std::sort(sorted_counts.begin(), sorted_counts.end(),
              [](const auto &a, const auto &b)
              {
                  return a.first > b.first;
              });

    std::cout << BLUE << "> " << RESET << "Fragment length counts:\n";
    for (const auto &[length, count] : sorted_counts)
    {
        std::cout << BLUE << "> " << RESET << "Length " << length << ": " << count << " occurrences\n";
        if (file_out)
        {
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
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solution is optimal.\n";
        // std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Objective value: " << response.objective_value() << std::endl;
    }
    else if (response.status() == operations_research::sat::CpSolverStatus::FEASIBLE)
    {
        std::cout << BLUE << "> " << RESET << "[" << elapsed_since_start() << "] Solution is not optimal, but feasible.\n";
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
