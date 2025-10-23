/**************************************************************/
/*         Now, from a solar system to another I fly          */
/*                I come with life in my veins                */
/*   Through darkened space I ran, I left a struggle behind   */
/*              And faced the fear that's alive               */
/*  I lost a world and traveled straight across from the sky  */
/*               To find the gates of the heart               */
/*                                                            */
/*                     To Sirius - Gojira                     */
/**************************************************************/

#include <tuple>
#include <queue>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "utils.hpp"
#include "solver.hpp"
#include "config.hpp"
#include "tables.hpp"
#include "logging.hpp"
#include "instance.hpp"
#include "interface.hpp"

#include "cxxopts.hpp"

#include "ortools/util/sigint.h"
#include "ortools/sat/cp_model.h"
#include "ortools/base/logging.h"
#include "ortools/base/init_google.h"

int main(int argc, char *argv[])
{
    // Systematische Identifikation Redundanter, Identisch Uebersetzter Sequenzen
    sirius::utils::preparse_quiet_flag(argc, argv);
    ::google::InitGoogleLogging("SIRIUS");
    sirius::logg.begruessung();  // Guten Tag

    operations_research::SigintHandler sigintHandler;
    sigintHandler.Register([]() {
        sirius::logg.on_sigint();
    });  // Pretty CTRL + C

    sirius::SIRIUSConfig config;
    sirius::SIRIUSTables tables;
    sirius::SIRIUSInstance instance;
    sirius::SIRIUSInterface interface;

    if (argc == 1)
    { 
        std::tie(instance, config) = interface.gather_inputs_interactively(tables); 
    }
    else
    { 
        std::tie(instance, config) = interface.gather_inputs_from_flags(argc, argv, tables); 
    }

    if (instance.hard_filter_by_rscu || instance.soft_filter_by_rscu)
    { 
        tables.build_rscu_map_from_csv(instance.codon_usage_path); 
    }

    sirius::SIRIUSSolver sirius_solver(instance, config, tables);

    // ----- GO! -----
    sirius_solver.init_new_model();
    sirius_solver.build_model();
    sirius_solver.solve_model();
    // ---------------

    // Variables for each base – Flatten sequence_vars_list
    std::vector<operations_research::sat::BoolVar> base_vars;
    for (const auto &a : sirius_solver.sequence_vars_list) {
        for (const auto &b : a) {
            for (const auto &c : b) {
                for (const auto &var : c) {
                    base_vars.push_back(var);
                }
            }
        }
    }

    std::queue<char> int_vars;
    for (const auto &var : base_vars)
    {
        if (operations_research::sat::SolutionIntegerValue(sirius_solver.response, var) == 1)
        {
            int_vars.push(var.Name().front());
        }
    }

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

    std::string output_folder = sirius::utils::create_output_folder(instance.output_folder);
    std::string summary_filename = output_folder + "/" + sirius::utils::timestamped_filename("summary");
    std::string sequences_filename = output_folder + "/" + sirius::utils::timestamped_filename("sequences");
    std::string length_counts_filename = output_folder + "/" + sirius::utils::timestamped_filename("length_counts");

    sirius::utils::validate_translated_proteins(all_out_seqs, instance.init_target_protein, tables);
    sirius::utils::write_sequences_to_file_and_console(all_out_seqs, sequences_filename);
    sirius::utils::write_summary_to_file(sirius::logg.summary, summary_filename);

    auto [
        all_stretches,
        length_counts
    ] = sirius::utils::find_all_homologous_stretches_and_count_lengths(all_out_seqs);

    try
    {
        std::ofstream out_lengths(length_counts_filename);
        sirius::utils::print_length_counts(length_counts, &out_lengths);
        out_lengths.close();
    }
    catch (const std::exception &e)
    {
        throw e;
    }

    return 0;
}