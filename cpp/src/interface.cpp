#include <utility>
#include <fstream>

#include "utils.hpp"
#include "config.hpp"
#include "tables.hpp"
#include "logging.hpp"
#include "instance.hpp"
#include "interface.hpp"

#include "cxxopts.hpp"

namespace sirius 
{
    std::pair<SIRIUSInstance, SIRIUSConfig>
    SIRIUSInterface::gather_inputs_from_flags(int argc, char* argv[], const SIRIUSTables& tables)
    {
        // ---- defaults ----
        int         num_sequences          = 2;
        std::string init_target_protein    = "MALEEINENSTERN";
        int         num_workers            = utils::available_threads_linux();
        bool        show_ortools_log       = false;
        double      relative_gap_limit     = 0.0;
        int         max_time_in_seconds    = 600;

        bool   user_set_alpha              = false;
        bool   user_set_rscu_ratio         = false;
        bool   user_set_gc_end_rscu_thresh = false;
        bool   user_set_hard_rscu_thresh   = false;
        bool   user_set_soft_rscu_thresh   = false;

        std::string codon_usage_path = "";
        // Let Python wrapper tell us the packaged default CSV (if available):
        std::string packaged_default_csv;
        if (const char* p = std::getenv("SIRIUS_DEFAULT_HOST_CSV"))
        {
            packaged_default_csv = p;
        } 
        else
        {
            packaged_default_csv = "";
        }

        std::string output_folder    = "sirius_out";

        double rscu_alpha            = 10.0;
        double max_low_rscu_ratio    = 0.3;
        double hard_rscu_threshold   = 0.0;
        double soft_rscu_threshold   = 0.0;
        double gc_end_rscu_threshold = 0.0;

        // ---- define CLI ----
        cxxopts::Options opts("sirius", "*** Flight Manual ***");

        opts.positional_help("<protein-sequence>");

        opts.add_options()
            // required
            ("p,prot", "Protein sequence (AA string). You may also provide it positionally.",
                cxxopts::value<std::string>()->default_value(init_target_protein))

            // knobs
            ("n,num", "Number of variants",
                cxxopts::value<int>()->default_value(std::to_string(num_sequences)))
            ("w,num-workers", "CPU workers (threads)",
                cxxopts::value<int>()->default_value(std::to_string(num_workers)))
            ("t,max-sec", "Time limit (seconds)",
                cxxopts::value<int>()->default_value(std::to_string(max_time_in_seconds)))
            ("R,relative-gap-limit", "Relative LIP gap limit (0..1)",
                cxxopts::value<double>()->default_value(std::to_string(relative_gap_limit)))
            ("v,show-ortools-log", "Enable solver logging",
                cxxopts::value<bool>()->default_value("false")->implicit_value("true"))

            ("o,output-folder", "Output folder name",
                cxxopts::value<std::string>()->default_value("sirius_out"))

            // RSCU / CSV
            ("c,codon-usage-fpath", "Path to codon usage CSV (overrides packaged default)",
                cxxopts::value<std::string>()->default_value(packaged_default_csv))
            ("h,hard-rscu-thresh", "Hard RSCU threshold",
                cxxopts::value<double>()->default_value("0.0"))
            ("s,soft-rscu-thresh", "Soft RSCU threshold",
                cxxopts::value<double>()->default_value("0.0"))
            // ("gc-end-rscu-thresh", "GC-end RSCU threshold",
            //     cxxopts::value<double>()->default_value("0.0"))
            ("a,rscu-alpha", "Alpha for soft RSCU",
                cxxopts::value<double>()->default_value("10.0"))
            ("m,max-low-rscu-ratio", "Max fraction (0..1) of low-RSCU codons allowed in soft filtering",
                cxxopts::value<double>()->default_value("0.3"))

            ("q,quiet", "Quiet mode")
            ("h,help", "Show help");

        // Make "prot" accepted positionally too
        opts.parse_positional({"prot"});

        cxxopts::ParseResult result;
        try 
        {
            result = opts.parse(argc, argv);
        }
        catch (const std::exception& e)
        {
            logg.throw_formatted_error(std::string("Error parsing arguments: ") +
                e.what() + "\n\n" + opts.help());
        }

        if (result.count("help"))
        {
            std::cout << opts.help() << std::endl;
            std::exit(0);
        }

        // ---- read values ----
        if (result.count("prot"))
        {
            init_target_protein = result["prot"].as<std::string>();

            if (init_target_protein[0] != 'M')
            {
                init_target_protein = "M" + init_target_protein;
                logg.print_warning_newline("Böse böse... Auto adding 'M' to start your protein :)");
            }
        } 
        else
        {
            logg.throw_formatted_error(
                "Error: missing protein sequence. Use -p/--prot or positional <protein-sequence>.\n\n" +
                opts.help());
        }

        num_sequences       = result["num"].as<int>();
        num_workers         = result["num-workers"].as<int>();
        max_time_in_seconds = result["max-sec"].as<int>();
        relative_gap_limit  = result["relative-gap-limit"].as<double>();
        show_ortools_log    = result["show-ortools-log"].as<bool>();

        codon_usage_path    = result["codon-usage-fpath"].as<std::string>();
        output_folder       = result["output-folder"].as<std::string>();

        hard_rscu_threshold    = result["hard-rscu-thresh"].as<double>();
        soft_rscu_threshold    = result["soft-rscu-thresh"].as<double>();
        // gc_end_rscu_threshold  = result["gc-end-rscu-thresh"].as<double>();
        rscu_alpha             = result["rscu-alpha"].as<double>();
        max_low_rscu_ratio     = result["max-low-rscu-ratio"].as<double>();

        user_set_hard_rscu_thresh   = result["hard-rscu-thresh"].count() && hard_rscu_threshold != 0.0;
        user_set_soft_rscu_thresh   = result["soft-rscu-thresh"].count() && soft_rscu_threshold != 0.0;
        // user_set_gc_end_rscu_thresh = result["gc-end-rscu-thresh"].count() && gc_end_rscu_threshold != 0.0;
        user_set_alpha              = result["rscu-alpha"].count() && rscu_alpha != 10.0;
        user_set_rscu_ratio         = result["max-low-rscu-ratio"].count() && max_low_rscu_ratio != 0.3;

        // ---- validations ----
        utils::validate_user_prot_input(init_target_protein, tables.reduced_codon_table);
        utils::validate_user_num_seq_input(num_sequences);
        utils::validate_csv(codon_usage_path, tables.invariant_codon_table);

        if (user_set_hard_rscu_thresh && user_set_soft_rscu_thresh) {
            logg.throw_formatted_error(
                "Error: Cannot both soft- and hard-filter RSCU. Define a value for only one.");
        }

        if ((user_set_hard_rscu_thresh || user_set_soft_rscu_thresh) && codon_usage_path.empty()) {
            logg.throw_formatted_error(
                "Error: Provide codon usage file (--codon-usage-fpath) when using RSCU filtering.");
        }

        if (user_set_soft_rscu_thresh && !user_set_rscu_ratio) {
            std::ostringstream oss;
            oss << "Warning: You did not specify an RSCU ratio for soft filtering. "
                << "Using default (" << std::fixed << std::setprecision(1) << max_low_rscu_ratio << ").";
            logg.print_warning_newline(oss.str());
        }

        if (user_set_soft_rscu_thresh && !user_set_alpha) {
            std::ostringstream oss;
            oss << "Warning: You did not specify an RSCU alpha for soft filtering. "
                << "Using default (" << std::fixed << std::setprecision(1) << rscu_alpha << ").";
            logg.print_warning_newline(oss.str());
        }

        if (!codon_usage_path.empty()) {
            utils::require_file_if_provided(codon_usage_path, "codon usage file");
        }

        // bounds checks
        if (max_low_rscu_ratio < 0.0 || max_low_rscu_ratio > 1.0) {
            logg.throw_formatted_error("Error: --max-low-rscu-ratio must be in [0,1].");
        }
        if (relative_gap_limit < 0.0 || relative_gap_limit > 1.0) {
            logg.throw_formatted_error("Error: --relative-gap-limit must be in [0,1].");
        }
        if (num_sequences <= 0) {
            logg.throw_formatted_error("Error: --num must be > 0.");
        }
        if (num_workers <= 0) {
            logg.throw_formatted_error("Error: --num-workers must be > 0.");
        }
        if (max_time_in_seconds <= 0) {
            logg.throw_formatted_error("Error: --max-sec must be > 0.");
        }

        utils::print_inputs(init_target_protein, num_sequences);

        SIRIUSInstance instance(
            num_sequences,
            init_target_protein,
            tables,
            codon_usage_path,         // may be empty -> GD uses packaged default via env
            output_folder,
            hard_rscu_threshold,
            soft_rscu_threshold,
            gc_end_rscu_threshold,
            rscu_alpha,
            max_low_rscu_ratio
        );

        SIRIUSConfig config(
            show_ortools_log,
            num_workers,
            /*linearization_level=*/3,
            relative_gap_limit,
            max_time_in_seconds
        );

        std::ostringstream sum;
        sum << "*** SIRIUS Summary ***\n"
            << "Protein:     " << init_target_protein << "\n"
            << "Variants:    " << num_sequences << "\n";

        if (user_set_hard_rscu_thresh)
        {
            sum << "Codon CSV:   " << codon_usage_path << "\n"
                << "Mode:        HARD RSCU\n"
                << "  hard RSCU: " << hard_rscu_threshold << "\n";
                // << "  GC-end:    " << gc_end_rscu_threshold << "\n";
        } 
        else if (user_set_soft_rscu_thresh) 
        {
            sum << "Codon CSV:   " << codon_usage_path << "\n"
                << "Mode:        SOFT RSCU\n"
                << "  soft RSCU: " << soft_rscu_threshold << "\n"
                << "  alpha:     " << rscu_alpha << "\n"
                << "  max ratio: " << max_low_rscu_ratio << "\n";
        }
        else 
        {
            sum << "Mode:        no RSCU filtering\n";
        }

        sum << "Timeout(s):  " << max_time_in_seconds << "\n"
            << "Threads:     " << num_workers;

        logg.summary = sum.str();

        return {instance, config};
    }

    std::pair<SIRIUSInstance, SIRIUSConfig>
    SIRIUSInterface::gather_inputs_interactively(const SIRIUSTables& tables)
    {
        // defaults
        const std::string default_prot  = "MALEEINENSTERN";
        const std::string output_folder = "sirius_out";
        int    num_sequences            = 2;
        int    max_time_in_seconds      = 600;
        int    num_workers              = utils::available_threads_linux();

        // RSCU knobs
        bool   use_hard_rscu            = false;
        bool   use_soft_rscu            = false;
        double hard_rscu_threshold      = 0.5;
        double soft_rscu_threshold      = 0.7;
        // double gc_end_rscu_threshold    = 0.5;
        double gc_end_rscu_threshold    = 0.0;
        double rscu_alpha               = 10.0;
        double max_low_rscu_ratio       = 0.3;

        // Let Python wrapper tell us the packaged default CSV (if available):
        std::string packaged_default_csv;
        if (const char* p = std::getenv("SIRIUS_DEFAULT_HOST_CSV"))
        {
            packaged_default_csv = p;
        } 
        else
        {
            packaged_default_csv = "";
        }

        // 1) Protein & count
        std::string init_target_protein = utils::ask_protein("Gather your protein: ", default_prot, tables);

        if (init_target_protein[0] != 'M')
        {
            init_target_protein = "M" + init_target_protein;
            logg.print_warning_newline("Böse böse... Auto adding 'M' to start your protein :)");
        }

        num_sequences = utils::ask_positive_int("And the number of sequences [Default 2]: ", 2);

        // 2) Choose filtering mode
        use_hard_rscu = utils::ask_yes_no("Hard filter by RSCU? [Default no]: ", false);
        if (!use_hard_rscu) 
        {
            use_soft_rscu = utils::ask_yes_no("Soft filter by RSCU? [Default no]: ", false);
        } 
        else 
        {
            use_soft_rscu = false;
        }

        std::string codon_usage_path; // empty -> use packaged default

        if (use_soft_rscu) 
        {
            soft_rscu_threshold    = utils::ask_number<double>("Soft RSCU threshold [Default 0.7]: ", 0.7);
            rscu_alpha             = utils::ask_number<double>("Soft RSCU alpha [Default 10.0]: ", 10.0);
            max_low_rscu_ratio     = utils::ask_number<double>("Max low-RSCU ratio (0..1) [Default 0.3]: ", 0.3, true, 0.0, 1.0);
            // CSV optional: prompt allows Enter to use packaged default
            codon_usage_path       = utils::ask_codon_csv_or_default(packaged_default_csv);
        }

        if (use_hard_rscu) 
        {
            hard_rscu_threshold     = utils::ask_number<double>("Hard RSCU threshold [Default 0.5]: ", 0.5);
            // gc_end_rscu_threshold   = utils::ask_number<double>("GC-ending RSCU threshold [Default 0.5]: ", 0.5);
            // CSV required for hard filtering; still allow Enter to use packaged default
            codon_usage_path        = utils::ask_codon_csv_or_default(packaged_default_csv);
        }

        if (use_hard_rscu || use_soft_rscu)
        {   
            if (codon_usage_path.empty())
            {
                codon_usage_path = packaged_default_csv;
            }
            else if (!codon_usage_path.empty())
            {
                if (!utils::file_exists(codon_usage_path))
                {
                    logg.throw_formatted_error("Error: codon usage file not found at '" + codon_usage_path + "'.");
                }
            }

            utils::validate_csv(codon_usage_path, tables.invariant_codon_table);
        }

        // 3) Solver settings (keep earlier fixed defaults, or ask interactively)
        // Uncomment to ask:
        // max_time_in_seconds = ask_positive_int("Max time in seconds [Default 600]: ", 600);
        // num_workers = ask_positive_int("CPU workers [Default auto]: ", available_threads_linux());

        // 4) Summary + confirm
        if (!logg.quiet)
        {
            std::ostringstream sum;
            sum << "*** SIRIUS Summary ***\n"
                << "Protein:     " << init_target_protein << "\n"
                << "Variants:    " << num_sequences << "\n";

            if (use_hard_rscu)
            {
                sum << "Codon CSV:   " << (codon_usage_path.empty() ? (packaged_default_csv + " (packaged)") : codon_usage_path) << "\n"
                    << "Mode:        HARD RSCU\n"
                    << "  hard RSCU: " << hard_rscu_threshold << "\n";
                    // << "  GC-end:    " << gc_end_rscu_threshold << "\n";
            } 
            else if (use_soft_rscu) 
            {
                sum << "Codon CSV:   " << (codon_usage_path.empty() ? (packaged_default_csv + " (packaged)") : codon_usage_path) << "\n"
                    << "Mode:        SOFT RSCU\n"
                    << "  soft RSCU: " << soft_rscu_threshold << "\n"
                    << "  alpha:     " << rscu_alpha << "\n"
                    << "  max ratio: " << max_low_rscu_ratio << "\n";
            }
            else 
            {
                sum << "Mode:        no RSCU filtering\n";
            }

            sum << "Timeout(s):  " << max_time_in_seconds << "\n"
                << "Threads:     " << num_workers;

            logg.print_summary(sum.str());

            if (!utils::ask_yes_no("Proceed? [Y/n]: ", true))
            {
                logg.throw_formatted_error("Aborted.");
            }
        }

        utils::print_inputs(init_target_protein, num_sequences);

        // Construct outputs
        SIRIUSInstance instance(
            num_sequences,
            init_target_protein,
            tables,
            codon_usage_path,
            output_folder,
            use_hard_rscu ? hard_rscu_threshold : 0.0,
            use_soft_rscu ? soft_rscu_threshold : 0.0,
            use_hard_rscu ? gc_end_rscu_threshold : 0.0,
            rscu_alpha,
            max_low_rscu_ratio
        );

        // Note: set 'show_ortools_log' to false here for interactive default
        SIRIUSConfig config(/*show_ortools_log=*/false,
                            /*num_workers=*/num_workers,
                            /*linearization_level=*/3,
                            /*relative_gap_limit=*/0.0,
                            /*max_time_in_seconds=*/max_time_in_seconds);

        return {instance, config};
    }

} // namespace sirius