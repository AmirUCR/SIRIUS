#include <vector>
#include <string>
#include <thread>
#include <cstdio>
#include <iomanip>
#include <filesystem>
#include <unordered_map>

#include "utils.hpp"
#include "logging.hpp"
#include "instance.hpp"

#include "csv.hpp"

namespace sirius::utils
{
    // Pre-scan argv for a boolean --quiet / -q before any printing/log setup.
    // Accepts: --quiet, --quiet=true/false, -q, -q=true/false
    void preparse_quiet_flag(int argc, char* argv[])
    {
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg(argv[i]);
            if (arg == "-q" || arg == "--quiet")
            {
                logg.quiet = true;
                return;
            }
            if (arg.rfind("--quiet=", 0) == 0)
            {
                const auto val = arg.substr(8);
                logg.quiet = (val == "1" || val == "true" || val == "TRUE" || val == "True");
                return;
            }
            if (arg.rfind("-q=", 0) == 0)
            {
                const auto val = arg.substr(3);
                logg.quiet = (val == "1" || val == "true" || val == "TRUE" || val == "True");
                return;
            }
        }
    }

    // Returns min(available cores, 32)
    // OR-Tools performs best with 32 threads
    // based on my experience.
    unsigned available_threads_linux() {
        unsigned n = 0;

        cpu_set_t set;
        CPU_ZERO(&set);
        if (sched_getaffinity(0, sizeof(set), &set) == 0) {
            unsigned m = CPU_COUNT(&set);
            if (m > 0) {
                n = m;
            }
        }

        if (n == 0) {
            n = std::thread::hardware_concurrency();
        }
        if (n == 0) {
            n = 1;
        }

        return std::min(n, 32u);
    }

    std::string ask_protein(
        const std::string& prompt,
        const std::string& def,
        const SIRIUSTables& tables)
    {
        while (true) 
        {
            std::string s;
            if (!prompt_line(prompt, s)) 
            {
                s.clear();
            }

            if (s.empty())
            {
                logg.print_info_newline("Setting protein to " + def);
                return def;
            }

            try 
            {
                validate_user_prot_input(s, tables.reduced_codon_table);
                return s;
            }
            catch (const std::exception& e) 
            {
                logg.print_warning_newline(e.what());
            }
        }
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
            logg.throw_formatted_error(oss.str());
        }
    }

    void validate_user_num_seq_input(int num_sequences)
    {
        // Check bounds on number of sequences
        if (num_sequences <= 1)
        {
            logg.throw_formatted_error("Must generate at least 2 sequences.");
        }
    }

    void validate_csv(
        const std::string &filename,
        const std::unordered_map<char, std::string> &invariant_codon_table)
    {
        csv::CSVReader reader(filename);

        double rscu;
        std::string codon;
        std::string aa_str;

        int n_rows = 0;
        // The 'reader' yields 'csv::CSVRow' objects.
        for (csv::CSVRow& row : reader)
        {
            try 
            {
                aa_str = row["AmOneLet"].get<std::string>();
            }
            catch (const std::exception& e) 
            {
                logg.throw_formatted_error("Error: Malformed data in " + 
                    filename + ", column AmOneLet, row " + std::to_string(n_rows + 1) + ". " + e.what());
            }

            try 
            {
                codon = row["Codon"].get<std::string>();
            }
            catch (const std::exception& e) 
            {
                logg.throw_formatted_error("Error: Malformed data in " + 
                    filename + ", column Codon, row " + std::to_string(n_rows + 1) + ". " + e.what());
            }

            try 
            {
                rscu = row["RSCU"].get<double>();
            }
            catch (const std::exception& e) 
            {
                logg.throw_formatted_error("Error: Malformed data in " + 
                    filename + ", column RSCU, row " + std::to_string(n_rows + 1) + ". " + e.what());
            }
            
            char aa = aa_str[0];
            auto it = invariant_codon_table.find(aa);
            if (it == invariant_codon_table.end() || aa_str.size() > 1)
            {
                logg.throw_formatted_error("Error: Unknown amino acid '" + aa_str + "'.");
            }

            if (codon.size() != 3)
            {
                logg.throw_formatted_error("Error: Malformed CSV line in " + 
                    filename + " (row " + std::to_string(n_rows) + "). Codon " + 
                    codon + " does not have three bases.");
            }

            if (aa_str.empty())
            {
                logg.throw_formatted_error("Error: Malformed CSV line in " + 
                    filename + " (row " + std::to_string(n_rows) + "). No amino acid found.");
            }

            ++n_rows;
        }
    }

    void validate_file_exists(const std::string &filename, const std::string &message)
    {
        if (!std::filesystem::exists(filename))
        {
            logg.throw_formatted_error("Error: file for " + message + " (\"" + filename + "\") does not exist.");
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

            logg.print_info_newline("Now, in that star's heart, forge " + std::to_string(num_sequences) + "x " + to_print);
        }
        else
        {
            logg.print_info_newline("Now, in that star's heart, forge " + std::to_string(num_sequences) + "x " + protein);
        }
    }

    // read a line; return false on EOF (so CI or piping won't hang)
    bool prompt_line(const std::string& prompt, std::string& out) 
    {
        logg.print_info(prompt);
        if (!std::getline(std::cin, out)) 
        {
            return false;
        }
        out = trim(out);
        return true;
    }

    // ask yes/no with default
    bool ask_yes_no(const std::string& prompt, bool default_val) {
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

            s = lower(s);
            if (s == "y" || s == "yes" || s == "true" || s == "1") 
            {
                return true;
            }
            if (s == "n" || s == "no"  || s == "false"|| s == "0")
            {
                return false;
            }
            logg.print_warning_newline("Please answer 'yes' or 'no'.");
        }
    }

    // ask int > 0
    int ask_positive_int(const std::string& prompt, int def) 
    {
        while (true) 
        {
            int v = ask_number<int>(prompt, def, true, 1, std::numeric_limits<int>::max());
            if (v > 0) 
            {
                return v;
            }
            logg.print_warning_newline("Value must be a positive integer.");
        }
    }

    std::string shell_quote(const std::string& s) {
        std::string out = "'";
        for (char c : s) {
            if (c == '\'') out += "'\"'\"'"; else out += c;
        }
        out += "'";
        return out;
    }

    std::string run_python(const std::string& command) {
        std::array<char, 512> buffer;
        std::string result;

        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            logg.throw_formatted_error("Dev Error: Failed to run Python command for GeneDiversifier.");
        }

        while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) 
        {
            result += buffer.data();
        }

        int rc = pclose(pipe);
        if (rc != 0) 
        {
            logg.throw_formatted_error(
                "Dev Error: GeneDiversifier failed to execute with code " + std::to_string(rc) + 
                ". Please report the steps to recreate this incident at https://github.com/ucrbioinfo/SIRIUS/issues.");
        }
        return result;
    }

    std::string call_gene_diversifier(
        const std::string& protein_seq,
        int n,
        std::string host_path,   // user can pass "", or a relative/absolute path
        double rscu_thresh,
        double gc_thresh)
    {
        const char* py_exec   = std::getenv("GENE_DIVERSIFIER_PYTHON");
        const char* py_script = std::getenv("GENE_DIVERSIFIER_PY");
        if (!py_exec || !py_script)
        {
            logg.throw_formatted_error("Dev Error: Environment variables for GeneDiversifier not set.");
        }

        // If user did not supply a path, use packaged default from the wrapper.
        if (host_path.empty())
        {
            const char* def_csv = std::getenv("SIRIUS_DEFAULT_HOST_CSV");

            if (!def_csv)
            {
                logg.throw_formatted_error("Dev Error: SIRIUS_DEFAULT_HOST_CSV not set.");
            }

            host_path = def_csv;
        }

        // Build safe command (quote each argument)
        std::ostringstream cmd;
        cmd << shell_quote(py_exec) << " "
            << shell_quote(py_script) << " "
            << shell_quote(protein_seq) << " "
            << shell_quote(std::to_string(n)) << " "
            << shell_quote(host_path) << " "
            << shell_quote(std::to_string(rscu_thresh)) << " "
            << shell_quote(std::to_string(gc_thresh))
            << " 2>&1"; // merge stderr to stdout for easier debugging

        std::string output = run_python(cmd.str());

        if (!output.empty() && output.back() == '\n')
        {
            output.pop_back();
        }
        return output;
    }

    // Helper: throw if a file path is non-empty but doesn't exist.
    void require_file_if_provided(const std::string& path, const char* what)
    {
        if (!path.empty())
        {
            std::ifstream f(path);
            if (!f.good())
            {
                std::ostringstream oss;
                oss << "Error: " << what << " not found at '" << path << "'.";
                logg.throw_formatted_error(oss.str());
            }
        }
    }

    // Returns either a validated user path or an empty string (meaning "use packaged default")
    std::string ask_codon_csv_or_default(const std::string& packaged_default_path)
    {
        std::ostringstream q;
        q << "Codon usage file path "
        << "[Enter for packaged default: " << packaged_default_path << "]: ";
        while (true) 
        {
            std::string s;
            if (!prompt_line(q.str(), s)) return ""; // non-interactive -> use default
            if (s.empty()) return "";                // use packaged default
            if (file_exists(s)) return s;
            logg.print_warning_newline("File not found: '" + s + "'. Try again or press Enter for default.");
        }
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
           << std::put_time(std::localtime(&now_c), "%Y-%m-%d_%H-%M-%S")
           << ".txt";
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

    void write_summary_to_file(
        const std::string &summary, 
        const std::string &base_filename)
    {
        std::string filename = generate_unique_filename(base_filename);
        std::ofstream out_file(filename);

        if (!out_file)
        {
            logg.throw_formatted_error("Error: Failed to open output file.");
        }

        out_file << summary << "\n";
        out_file.close();
    }

    void write_sequences_to_file_and_console(
        const std::vector<std::string> &sequences,
        const std::string &base_filename)
    {
        std::string filename = generate_unique_filename(base_filename);
        std::ofstream out_file(filename);

        if (!out_file)
        {
            logg.throw_formatted_error("Error: Failed to open output file.");
        }

        bool too_large = false;
        if (sequences.at(0).size() > 80)
        {
            too_large = true;
            logg.print_info_newline("Resulting sequences too long to display here.");
        }

        for (const auto &seq : sequences)
        {
            if (!too_large)
            {
                logg.print_info_newline(seq);
            }

            out_file << seq << "\n";
        }

        out_file.close();
        logg.print_info_newline("Find your sequences in " + filename);
    }

    // Function to validate protein translation of sequences
    void validate_translated_proteins(
        const std::vector<std::string> &sequences,
        const std::string &target_protein,
        const SIRIUSTables &tables)
    {
        for (size_t i = 0; i < sequences.size(); ++i)
        {
            std::string protein = translate_dna_to_protein(sequences[i], tables);

            if (protein != target_protein)
            {
                logg.print_error_newline("Dev error: Sequence " + std::to_string(i) + " and " + sequences[i] + " are not the same!");
            }
        }
    }

    std::string translate_dna_to_protein(
        const std::string &dna,
        const SIRIUSTables &tables)
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

    void print_length_counts(
        const std::unordered_map<int, int> &length_counts, 
        std::ostream *file_out)
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

        logg.print_info_newline("Fragment length counts:");
        for (const auto &[length, count] : sorted_counts)
        {
            logg.print_info_newline("Length " + std::to_string(length) + ": " + std::to_string(count) + " occurrences");
            if (file_out)
            {
                *file_out << "Length " << length << ": " << count << " occurrences\n";
            }
        }
    }

    // Print out solver response message
    void check_response(const operations_research::sat::CpSolverResponse &response)
    {
        if (response.status() == operations_research::sat::CpSolverStatus::OPTIMAL)
        {
            logg.print_timed_info_newline("Solution is optimal.");
            logg.print_timed_info_newline("Objective value: " + std::to_string(response.objective_value()));
        }
        else if (response.status() == operations_research::sat::CpSolverStatus::FEASIBLE)
        {
            logg.print_timed_info_newline("Solution is not proven optimal, but feasible.");
            logg.print_timed_info_newline("Objective value: " + std::to_string(response.objective_value()));
        }
        else if (response.status() == operations_research::sat::CpSolverStatus::INFEASIBLE)
        {
            logg.print_timed_error_newline("DEV ERROR: Infeasible. Please report the steps to recreate this incident at https://github.com/ucrbioinfo/SIRIUS/issues.");
        }
        else if (response.status() == operations_research::sat::CpSolverStatus::MODEL_INVALID)
        {
            logg.print_timed_error_newline("DEV ERROR: MODEL INVALID. Please report the steps to recreate this incident at https://github.com/ucrbioinfo/SIRIUS/issues.");
        }
    }
}
