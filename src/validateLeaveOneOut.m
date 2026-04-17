function loo_results = validateLeaveOneOut(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, ...
        TF_list, z0_base, TF_exp, beta_matrix, gene_to_rxn_map, tau, ...
        resultsTableFilter, final_results, params)

    fprintf('\n--- Starting leave-one-out validation ---\n');

    loo_results = table();

    if isempty(resultsTableFilter) || isempty(final_results)
        return;
    end

    [nMets_loo, nRxns_loo] = size(modelWT_FBA.S);
    v0_loo = 0;
    a0_loo = nRxns_loo;
    b0_loo = 2 * nRxns_loo;
    nVars_loo = 3 * nRxns_loo;
    vtype_loo = repmat('C', nVars_loo, 1);

    if ~isfield(params, 'flexBiomass') || isempty(params.flexBiomass)
        params.flexBiomass = 0.5;
    end

    gparams = struct();
    gparams.OutputFlag = 0;

    for s = 1:height(resultsTableFilter)
        if ismember('Strategy', resultsTableFilter.Properties.VariableNames)
            strat_name = string(resultsTableFilter.Strategy{s});
        else
            strat_name = "Multi_Iter" + string(resultsTableFilter.Iteration(s));
        end

        z_multi = resultsTableFilter.z_sol{s};
        kappa_val = resultsTableFilter.Kappa(s);

        if ismember('Max_product', final_results.Properties.VariableNames)
            orig_prod = final_results.Max_product(s);
        elseif ismember('Product', final_results.Properties.VariableNames)
            orig_prod = final_results.Product(s);
        else
            error('final_results must contain either "Max_product" or "Product".');
        end

        ko_indices  = find(z0_base == 1 & z_multi < 0.5);
        act_indices = find(z0_base == 0 & z_multi > 0.5);
        all_perturb_indices = [ko_indices; act_indices];

        if isempty(all_perturb_indices)
            continue;
        end

        fprintf('Analyzing %s (%d total perturbations)...\n', strat_name, length(all_perturb_indices));

        for p = 1:length(all_perturb_indices)
            tf_idx = all_perturb_indices(p);
            tf_name = string(TF_list{tf_idx});

            % Leave one out by restoring this TF to its WT state
            z_loo = z_multi;
            z_loo(tf_idx) = z0_base(tf_idx);

            perturb_type = "KO_Restored";
            if z0_base(tf_idx) == 0
                perturb_type = "Activation_Restored";
            end

            TF_activity_loo = TF_exp .* z_loo;
            y_loo = beta_matrix * TF_activity_loo;
            y_loo = gene_to_rxn_map * y_loo;
            y_loo = y_loo ./ tau;
            y_loo(isnan(y_loo)) = 1;
            y_loo(~isfinite(y_loo)) = 1;

            A_loo = sparse(0, nVars_loo);
            rhs_loo = [];
            sen_loo = '';

            A_loo = [A_loo; sparse(modelWT_FBA.S), sparse(nMets_loo, 2*nRxns_loo)];
            rhs_loo = [rhs_loo; zeros(nMets_loo,1)];
            sen_loo = [sen_loo; repmat('=', nMets_loo, 1)];

            for i = 1:nRxns_loo
                y_i = y_loo(i);
                if any(gene_to_rxn_map(i,:)) && y_i < 0.99
                    row_lb = zeros(1, nVars_loo);
                    row_lb(v0_loo + i) = 1;
                    row_lb(a0_loo + i) = 1;
                    A_loo = [A_loo; row_lb];
                    rhs_loo = [rhs_loo; modelWT_FBA.lb(i) * y_i];
                    sen_loo = [sen_loo; '>'];

                    row_ub = zeros(1, nVars_loo);
                    row_ub(v0_loo + i) = 1;
                    row_ub(b0_loo + i) = -1;
                    A_loo = [A_loo; row_ub];
                    rhs_loo = [rhs_loo; modelWT_FBA.ub(i) * y_i];
                    sen_loo = [sen_loo; '<'];
                end
            end

            lb_loo = [modelWT_FBA.lb; zeros(2*nRxns_loo,1)];
            ub_loo = [modelWT_FBA.ub; inf(2*nRxns_loo,1)];
            lb_loo(v0_loo + biomass_idx) = biomassWT_FBA * params.flexBiomass;
            ub_loo(v0_loo + biomass_idx) = 1000;

            model_loo = struct('A', A_loo, 'rhs', rhs_loo, 'sense', sen_loo, ...
                               'lb', lb_loo, 'ub', ub_loo, 'vtype', vtype_loo, 'modelsense', 'max');

            % Reward product heavily to overcome kappa slack penalties
            weight_prod = max(1e4, kappa_val * 10);
            model_loo.obj = [zeros(nRxns_loo,1); -kappa_val*ones(nRxns_loo,1); -kappa_val*ones(nRxns_loo,1)];
            model_loo.obj(v0_loo + prod_idx) = weight_prod;

            result_loo = gurobi(model_loo, gparams);

            if isfield(result_loo, 'status') && strcmp(result_loo.status, 'OPTIMAL')
                v_sol_loo = result_loo.x(1:nRxns_loo);
                loo_prod = v_sol_loo(prod_idx);
            else
                loo_prod = NaN;
            end

            if ~isnan(loo_prod) && orig_prod > 1e-6
                prod_drop_pct = 100 * (orig_prod - loo_prod) / orig_prod;
            else
                prod_drop_pct = 0;
            end

            loo_results = [loo_results; ...
                table({char(strat_name)}, {char(tf_name)}, {char(perturb_type)}, orig_prod, loo_prod, prod_drop_pct, ...
                'VariableNames', {'Strategy', 'Restored_TF', 'Action', 'Multi_Product', 'LOO_Product', 'Product_Drop_Pct'})]; %#ok<AGROW>

            fprintf('  Restored %s: Prod went from %.4f -> %.4f (Drop: %5.1f%%) \n', ...
                tf_name, orig_prod, loo_prod, prod_drop_pct);
        end
    end
end
