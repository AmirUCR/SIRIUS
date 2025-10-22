#pragma once
#include <map>
#include <vector>
#include <fstream>
#include <utility>
#include <cstring>
#include <algorithm>
#include <unordered_map>

#include "config.hpp"
#include "tables.hpp"
#include "logging.hpp"
#include "instance.hpp"

#include "ortools/sat/cp_model_solver.h"

// Type alias for a homology stretch (start index, end index, length)
using Stretch = std::tuple<int, int, int>;

namespace sirius::utils
{
    void preparse_quiet_flag(int argc, char* argv[]);

    inline bool str_ieq(const char* a, const char* b)
    {
        return std::strcmp(a, b) == 0;
    }

    unsigned available_threads_linux();

    std::string shell_quote(const std::string& s);
    
    std::string run_python(const std::string& command);

    std::string ask_protein(
        const std::string& prompt,
        const std::string& def,
        const SIRIUSTables& tables);

    bool prompt_line(const std::string& prompt, std::string& out);

    template <typename T>
    T ask_number(
        const std::string& prompt,
        T default_val,
        bool enforce_range=false, 
        T minv=T{},
        T maxv=T{})
    {
        while (true)
        {
            std::string s;
            if (!prompt_line(prompt, s))
            {
                return default_val;
            }
            if (s.empty()) 
            {
                return default_val;
            }

            std::istringstream iss(s);
            T v;
            iss >> v;

            if (!iss.fail() && iss.eof()) 
            {
                if (!enforce_range || (v >= minv && v <= maxv))
                {
                    return v;
                }
                std::ostringstream oss;
                oss << "Value must be in [" << minv << ", " << maxv << "].";
                logg.print_warning_newline(oss.str());
            } 
            else 
            {
                logg.print_warning_newline("Please enter a valid number.");
            }
        }
    }

    void validate_user_prot_input(
        const std::string &protein,
        const std::unordered_map<char, std::vector<std::string>> &reduced_codon_table);

    void validate_user_num_seq_input(int num_sequences);

    void validate_csv(
        const std::string &filename, 
        const std::unordered_map<char, std::string> &invariant_codon_table);

    void validate_file_exists(const std::string &filename, const std::string &message);

    void print_inputs(const std::string &protein, int num_sequences);

    inline std::string trim(std::string s) 
    {
        auto notspace = [](int ch){ return !std::isspace(ch); };
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), notspace));
        s.erase(std::find_if(s.rbegin(), s.rend(), notspace).base(), s.end());
        return s;
    }

    inline std::string lower(std::string s) 
    {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        return s;
    }

    inline bool file_exists(const std::string& path) 
    {
        std::ifstream f(path);
        return f.good();
    }

    bool ask_yes_no(const std::string& prompt, bool default_val);

    int ask_positive_int(const std::string& prompt, int def);

    std::string call_gene_diversifier(
        const std::string& protein_seq,
        int n,
        std::string host_path,   // user can pass "", or a relative/absolute path
        double rscu_thresh,
        double gc_thresh);

    std::pair<SIRIUSInstance, SIRIUSConfig>
    gather_inputs_from_flags(int argc, char* argv[], const SIRIUSTables& tables);

    std::pair<SIRIUSInstance, SIRIUSConfig>
    gather_inputs_interactively(SIRIUSTables& tables);
    
    std::pair<std::map<std::string, std::vector<Stretch>>, std::unordered_map<int, int>>
    find_all_homologous_stretches_and_count_lengths(const std::vector<std::string> &sequences);

    std::vector<Stretch> find_homologous_stretches(const std::string &seq1, const std::string &seq2);
    
    void require_file_if_provided(const std::string& path, const char* what);

    std::string ask_codon_csv_or_default(const std::string& packaged_default_path);

    std::string elapsed_since_start();

    std::string create_output_folder(const std::string &base);

    void write_summary_to_file(
        const std::string &summary, 
        const std::string &base_filename);

    std::string timestamped_filename(const std::string &prefix);

    std::string generate_unique_filename(const std::string &base_name);

    void write_sequences_to_file_and_console(
        const std::vector<std::string> &sequences,
        const std::string &base_filename);

    void validate_translated_proteins(
        const std::vector<std::string> &sequences,
        const std::string &target_protein,
        const SIRIUSTables &tables);

    std::string translate_dna_to_protein(
        const std::string &dna,
        const SIRIUSTables &tables);

    void print_length_counts(
        const std::unordered_map<int, int> &length_counts, 
        std::ostream *file_out);

    void check_response(const operations_research::sat::CpSolverResponse &response);
} // namespace sirius
