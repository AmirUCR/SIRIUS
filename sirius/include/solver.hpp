#pragma once
#include <string>
#include <vector>
#include <unordered_map>

#include "config.hpp"
#include "tables.hpp"
#include "instance.hpp"

#include "ortools/sat/cp_model.h"

namespace sirius
{
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

        int dna_size = 0;
        int codon_len = 3;
        int decidable_protein_length = 0;

        std::string dna_with_holes = "";
        std::string decidable_protein = "";

        std::vector<operations_research::sat::BoolVar> all_vars;
        std::vector<std::vector<operations_research::sat::LinearExpr>> all_pairs_y_terms;
        std::vector<std::vector<operations_research::sat::BoolVar>> all_pairs_z_terms;
        std::vector<operations_research::sat::IntVar> objective;
        std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> sequence_vars_list;
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> all_obj_mults;
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;
        std::unordered_map<std::string, int> map_var_name_to_val;

        SIRIUSSolver(
            SIRIUSInstance instance, 
            SIRIUSConfig config, 
            SIRIUSTables tables);

        void assign_var_values_from_solution();
        void add_hints();
        void init_new_model();
        void build_model();

        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>
            create_obj_mults_vector(int max_frag_size);

        void set_minimize_objective_value(int max_frag_size);
        void solve_model();

        std::vector<std::vector<std::string>> generate_sequence_codons_with_rscu();
        std::vector<std::vector<std::string>> generate_sequence_codons_with_soft_rscu();
        std::vector<std::vector<std::string>> generate_sequence_codons_with_hard_rscu();

        std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>>
            add_base_variables_from_rscu();

        std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>>
            add_base_variables_from_prot();

        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>
            add_codon_constraints();

        std::vector<std::vector<operations_research::sat::BoolVar>>
            create_z_chained_vars();

        std::vector<std::vector<operations_research::sat::BoolVar>>
            generate_combinations(
                const std::vector<operations_research::sat::BoolVar> &vars,
                const int combination_size);

        void set_larger_fragment_obj_vals_to_zero(int max_frag_size);
        void create_additive_objective(int max_frag_size);
    };
}
