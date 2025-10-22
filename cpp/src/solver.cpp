#include <random>
#include <string>
#include <vector>
#include <iomanip>
#include <unordered_map>

#include "utils.hpp"
#include "config.hpp"
#include "tables.hpp"
#include "solver.hpp"
#include "logging.hpp"
#include "instance.hpp"

#include "ortools/sat/cp_model.h"

namespace sirius 
{
    SIRIUSSolver::SIRIUSSolver(
        SIRIUSInstance instance, 
        SIRIUSConfig config, 
        SIRIUSTables tables)
    : instance(std::move(instance))
    , config(std::move(config))
    , tables(std::move(tables))
    , model(std::make_unique<operations_research::sat::Model>())
    , cp_model()
    , gen(std::random_device{}()) // gen(42), // todo make injectable
    , dis(0.0, 1.0)
    {
        this->codon_len = 3;
        this->dna_size = this->instance.init_target_protein.size() * 3;
        this->dna_with_holes = this->instance.dna_with_holes;
        this->decidable_protein = this->instance.decidable_protein;
        this->decidable_protein_length = this->instance.decidable_protein_length;
    }

    void SIRIUSSolver::assign_var_values_from_solution()
    {
        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                int colcount = 0;
                std::vector<std::string> this_pair_z_pairs;

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
                        char base_x = codons_vec.at(codon_x_i).at(base_idx);
                        char base_y = codons_vec.at(codon_y_i).at(base_idx);

                        std::string var_name_x = absl::StrFormat(
                            "%c%d%d%d%d",
                            base_x,
                            s,
                            aa_pos_i,
                            codon_x_i,
                            base_idx);

                        std::string var_name_y = absl::StrFormat(
                            "%c%d%d%d%d",
                            base_y,
                            t,
                            aa_pos_i,
                            codon_y_i,
                            base_idx);

                        this->map_var_name_to_val[var_name_x] = 1;
                        this->map_var_name_to_val[var_name_y] = 1;

                        if (base_x == base_y)
                        {
                            std::string combo_name = var_name_x + var_name_y;

                            this->map_var_name_to_val[combo_name] = 1;
                            this_pair_z_pairs.push_back(combo_name);
                        }
                    }
                }

                std::string var_name = absl::StrFormat("z%d_p%d:%d", colcount, s, t);
                ++colcount;
                this->map_var_name_to_val[var_name] = 1;
            }
        }
    }

    void SIRIUSSolver::add_hints()
    {
        for (const operations_research::sat::BoolVar &v_old : all_vars)
        {
            auto it = map_var_name_to_val.find(v_old.Name());
            if (it == map_var_name_to_val.end())
            {
                continue;
            }

            operations_research::sat::BoolVar v_clone = cp_model.GetBoolVarFromProtoIndex(v_old.index());
            this->cp_model.AddHint(v_clone, it->second);
        }
    }

    void SIRIUSSolver::init_new_model()
    {
        cp_model = operations_research::sat::CpModelBuilder();
        model = std::make_unique<operations_research::sat::Model>();
        model->Add(operations_research::sat::NewSatParameters(this->config.parameters));
    }

    void SIRIUSSolver::build_model()
    {
        logg.print_timed_info_newline("Construct a model...");
        int max_frag_size = this->instance.warm_start_largest_fragment_length;

        this->sequence_vars_list =
            (this->instance.hard_filter_by_rscu || this->instance.soft_filter_by_rscu)
                ? add_base_variables_from_rscu()
                : add_base_variables_from_prot();
                
        this->sequence_codons_list = add_codon_constraints();
        this->all_pairs_z_terms = create_z_chained_vars();
        this->all_obj_mults = create_obj_mults_vector(max_frag_size);
        assign_var_values_from_solution();
        add_hints();
        set_larger_fragment_obj_vals_to_zero(max_frag_size);
        create_additive_objective(max_frag_size);
        set_minimize_objective_value(max_frag_size);
    }

    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>
        SIRIUSSolver::create_obj_mults_vector(int max_frag_size)
    {
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> all_obj_mults;
        all_obj_mults.resize(max_frag_size + 1);

        return all_obj_mults;
    }

    void SIRIUSSolver::set_minimize_objective_value(int max_frag_size)
    {
        std::vector<operations_research::sat::BoolVar> sum_vec;
        for (int i = max_frag_size - 1; i >= 0; --i)
        {
            for (const std::vector<operations_research::sat::BoolVar> &vec : this->all_obj_mults[i])
            {
                for (const operations_research::sat::BoolVar &var : vec)
                {
                    sum_vec.push_back(var);
                }
            }
        }

        this->cp_model.Minimize(operations_research::sat::LinearExpr::Sum(sum_vec));
    }

    void SIRIUSSolver::solve_model()
    {
        logg.print_timed_info_newline("Solve this enigma...");
        this->response = SolveCpModel(cp_model.Build(), model.get());
        utils::check_response(this->response);
    }

    std::vector<std::vector<std::string>>
    SIRIUSSolver::generate_sequence_codons_with_soft_rscu()
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

    std::vector<std::vector<std::string>> SIRIUSSolver::generate_sequence_codons_with_hard_rscu()
    {
        std::vector<std::vector<std::string>> sequence_codons;

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

                double rscu = this->tables.rscu_map.at(aa).at(full_codon);

                if (rscu > largest_rscu)
                {
                    largest_rscu = rscu;
                    codon_w_largest_rscu = full_codon;
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

                codons.push_back(full_codon);
                at_least_one_passed = true;
            }

            // If nothing passes the threshold, add the codon w/ highest RSCU
            if (!at_least_one_passed)
            {
                std::ostringstream msg;
                msg << std::fixed << std::setprecision(2);
                msg << "Warning: None of the codons for amino acid '" << aa
                    << "' (position " << i << ", starting from 0) pass RSCU threshold "
                    << this->instance.hard_rscu_threshold
                    << ". Using the codon with the highest RSCU ("
                    << codon_w_largest_rscu << ", " << largest_rscu << ").";

                logg.print_warning_newline(msg.str());

                codons.push_back(codon_w_largest_rscu);
            }

            sequence_codons.push_back(codons);
        }

        return sequence_codons;
    }

    std::vector<std::vector<std::string>> SIRIUSSolver::generate_sequence_codons_with_rscu()
    {
        return (this->instance.soft_filter_by_rscu) ? generate_sequence_codons_with_soft_rscu() : generate_sequence_codons_with_hard_rscu();
    }

    std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> 
        SIRIUSSolver::add_base_variables_from_rscu()
    {
        // Variables for each base
        std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> sequence_vars_list;
        std::vector<std::vector<std::string>> sequence_with_codons = generate_sequence_codons_with_rscu();

        for (int sequence_n = 0; sequence_n < this->instance.n; ++sequence_n)
        {
            std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> this_sequence_vars_list_of_list;

            for (int amino_acid_position = 0; amino_acid_position < sequence_with_codons.size(); ++amino_acid_position)
            {
                std::vector<std::vector<operations_research::sat::BoolVar>> codon_vars_list;

                for (size_t codon_number = 0; codon_number < sequence_with_codons.at(amino_acid_position).size(); ++codon_number)
                {
                    const std::string &codon = sequence_with_codons.at(amino_acid_position).at(codon_number);

                    std::vector<operations_research::sat::BoolVar> base_vars_list;
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

    std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> 
        SIRIUSSolver::add_base_variables_from_prot()
    {
        // Variables for each base
        std::vector<std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>>> sequence_vars_list;

        for (int sequence_n = 0; sequence_n < this->instance.n; ++sequence_n)
        {
            std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> this_sequence_vars_list_of_list;

            for (int amino_acid_position = 0; amino_acid_position < this->decidable_protein_length; ++amino_acid_position)
            {
                char amino_acid = this->decidable_protein[amino_acid_position];
                std::vector<std::vector<operations_research::sat::BoolVar>> codon_vars_list;

                for (size_t codon_number = 0; codon_number < this->tables.reduced_codon_table.at(amino_acid).size(); ++codon_number)
                {
                    const std::string &codon = this->tables.reduced_codon_table.at(amino_acid).at(codon_number);

                    std::vector<operations_research::sat::BoolVar> base_vars_list;
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

    std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> 
        SIRIUSSolver::add_codon_constraints()
    {
        std::vector<std::vector<std::vector<operations_research::sat::BoolVar>>> sequence_codons_list;

        // Constraints for Valid Codons
        for (const auto &this_sequence_vars_list_of_list_it : this->sequence_vars_list)
        {
            std::vector<std::vector<operations_research::sat::BoolVar>> this_sequence_codons;

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

                    operations_research::sat::BoolVar mult = this->cp_model.NewBoolVar().WithName(var_name);

                    const int group_size = codon_vars_list.size();
                    operations_research::sat::LinearExpr group_sum = operations_research::sat::LinearExpr::Sum(codon_vars_list);
                    // Enforce bi-directional implication
                    this->cp_model.AddEquality(group_sum, group_size).OnlyEnforceIf(mult);
                    this->cp_model.AddLessThan(group_sum, group_size).OnlyEnforceIf(mult.Not());

                    this->cp_model.AddEquality(operations_research::sat::LinearExpr::Sum(codon_vars_list), codon_vars_list.size() * mult);

                    codon_mult_list.push_back(mult);
                    this->all_vars.push_back(mult);
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

    std::vector<std::vector<operations_research::sat::BoolVar>> 
        SIRIUSSolver::create_z_chained_vars()
    {
        std::vector<std::vector<operations_research::sat::BoolVar>> all_pairs_z_terms;

        for (int s = 0; s < this->instance.n; ++s)
        {
            for (int t = s + 1; t < this->instance.n; ++t)
            {
                int colcount = 0;
                std::vector<operations_research::sat::BoolVar> this_pair_z_terms;

                for (int aa_pos_i = 0; aa_pos_i < this->decidable_protein_length; ++aa_pos_i)
                {
                    for (int codon_i = 0; codon_i < this->tables.reduced_codon_lengths_table.at(this->decidable_protein.at(aa_pos_i)); ++codon_i)
                    {
                        // int codon_positions = this->tables.reduced_codon_table.at(this->decidable_protein.at(aa_pos_i)).size();
                        int codon_positions = this->sequence_vars_list.at(s).at(aa_pos_i).size();
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

                                    this->cp_model.AddMultiplicationEquality(z, x, y);

                                    // Add z to z_terms
                                    z_terms += z;
                                    this->all_vars.push_back(z);
                                }
                            }
                        }

                        std::string var_name = absl::StrFormat("z%d_p%d:%d", colcount, s, t);
                        operations_research::sat::BoolVar z_j = cp_model.NewBoolVar().WithName(var_name);
                        colcount += 1;

                        this->cp_model.AddEquality(z_j, z_terms);

                        this_pair_z_terms.push_back(z_j);
                        this->all_vars.push_back(z_j);
                    }
                }
                all_pairs_z_terms.push_back(this_pair_z_terms);
            }
        }

        return all_pairs_z_terms;
    }

    // Create sliding window combinations out of vars depending on combination size
    // E.g., Input: vars = [z1, z2, z3], combination_size = 2
    //       Output: [[z1, z2], [z2, z3]]
    std::vector<std::vector<operations_research::sat::BoolVar>> 
        SIRIUSSolver::generate_combinations(
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

    void SIRIUSSolver::set_larger_fragment_obj_vals_to_zero(int max_frag_size)
    {
        if (max_frag_size >= this->decidable_protein_length)
            return;

        for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
        {
            std::vector<operations_research::sat::BoolVar> objective_vars;

            const std::vector<operations_research::sat::BoolVar> &pair_z_terms = this->all_pairs_z_terms.at(i);
            std::vector<std::vector<operations_research::sat::BoolVar>> z_combos = generate_combinations(
                pair_z_terms, max_frag_size + 1);

            for (const std::vector<operations_research::sat::BoolVar> &vi : z_combos)
            {
                std::string var_name = "";

                for (const operations_research::sat::BoolVar &i : vi)
                {
                    var_name += i.Name() + ".";
                }

                operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName("mult_" + var_name);

                // ----------------------------
                // mult ⇒ vi  (for all i)
                for (auto &v : vi)
                {
                    this->cp_model.AddImplication(mult, v);
                }
                // (AND(vi)) ⇒ mult  encoded as a single clause: (¬v1 ∨ ¬v2 ∨ ... ∨ mult)
                std::vector<operations_research::sat::BoolVar> clause;
                clause.reserve(vi.size() + 1);
                clause.push_back(mult);
                for (auto &v : vi)
                {
                    clause.push_back(v.Not());
                }
                this->cp_model.AddBoolOr(clause);
                // -----------------------------

                // this->all_vars.push_back(mult);
                objective_vars.push_back(mult);
            }

            this->all_obj_mults[max_frag_size].push_back(objective_vars);
            this->cp_model.AddEquality(
                operations_research::sat::LinearExpr::Sum(this->all_obj_mults[max_frag_size].back()), 0
            );
        }
    }

    void SIRIUSSolver::create_additive_objective(int max_frag_size)
    {
        for (int fragment_len = max_frag_size; fragment_len > 0; --fragment_len)
        {
            for (int i = 0; i < this->all_pairs_z_terms.size(); ++i)
            {
                std::vector<operations_research::sat::BoolVar> objective_vars;

                const std::vector<operations_research::sat::BoolVar> &pair_z_terms = this->all_pairs_z_terms.at(i);
                std::vector<std::vector<operations_research::sat::BoolVar>> z_combos = generate_combinations(
                    pair_z_terms,
                    fragment_len
                );

                for (const std::vector<operations_research::sat::BoolVar> &vi : z_combos)
                {
                    std::string var_name = "";

                    for (const operations_research::sat::BoolVar &i : vi)
                    {
                        var_name += i.Name() + ".";
                    }

                    operations_research::sat::BoolVar mult = cp_model.NewBoolVar().WithName("mult_" + var_name);

                    // ----------------------------
                    // mult ⇒ vi  (for all i)
                    for (auto &v : vi)
                    {
                        this->cp_model.AddImplication(mult, v);
                    }
                    // (AND(vi)) ⇒ mult  encoded as a single clause: (¬v1 ∨ ¬v2 ∨ ... ∨ mult)
                    std::vector<operations_research::sat::BoolVar> clause;
                    clause.reserve(vi.size() + 1);
                    clause.push_back(mult);
                    for (auto &v : vi)
                    {
                        clause.push_back(v.Not());
                    }
                    this->cp_model.AddBoolOr(clause);
                    // -----------------------------

                    // this->all_vars.push_back(mult);
                    objective_vars.push_back(mult);
                }

                this->all_obj_mults[fragment_len - 1].push_back(objective_vars);
            }
        }
    }

}   