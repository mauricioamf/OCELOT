function sKO_results = validateSingleKOs(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, ...
        TF_list, z0_base, TF_exp, beta_matrix, gene_to_rxn_map, tau, ...
        resultsTableFilter, final_results, params)

    if nargin < 13
        error('validateSingleKOs requires TF_exp, beta_matrix, gene_to_rxn_map, tau, resultsTableFilter, and final_results.');
    end

    if nargin < 14 || isempty(params)
        params = struct();
    end

    flexBiomass = getParam(params, 'flexBiomass', 0.5);
    nKOsLimit = getParam(params, 'nKOsLimit', 20);
    weightMin = getParam(params, 'singleKOWeightMin', 1e4);
    weightMultiplier = getParam(params, 'singleKOWeightMultiplier', 10);
    classPct = getParam(params, 'singleKOClassificationPct', 5);

    fprintf('\n--- Starting single TF KO Analysis ---\n');

    emptyVarTypes = {'cell','cell','cell','double','double','double'};
    emptyVarNames = {'Strategy','Perturbed_TF','Action','Multi_Product', ...
                     'sKO_Product','Product_Drop_Pct'};

    if isempty(resultsTableFilter) || height(resultsTableFilter) == 0 || ...
       isempty(final_results) || height(final_results) == 0
        sKO_results = table('Size', [0 numel(emptyVarNames)], ...
            'VariableTypes', emptyVarTypes, 'VariableNames', emptyVarNames);
        return;
    end

    [nMets_sKO, nRxns_sKO] = size(modelWT_FBA.S);
    v0_sKO = 0;
    a0_sKO = nRxns_sKO;
    b0_sKO = 2 * nRxns_sKO;
    nVars_sKO = 3 * nRxns_sKO;
    vtype_sKO = repmat('C', nVars_sKO, 1);

    TF_list = toCellArrayLocal(TF_list);

    sKO_results = table('Size', [0 numel(emptyVarNames)], ...
        'VariableTypes', emptyVarTypes, 'VariableNames', emptyVarNames);

    for s = 1:height(resultsTableFilter)
        strat_name = sprintf('Multi_Iter%d', resultsTableFilter.Iteration(s));
        z_multi = resultsTableFilter.z_sol{s};
        kappa_val = resultsTableFilter.Kappa(s);

        if s <= height(final_results) && ismember('Max_product', final_results.Properties.VariableNames)
            orig_prod = final_results.Max_product(s);
        elseif ismember('Product', resultsTableFilter.Properties.VariableNames)
            orig_prod = resultsTableFilter.Product(s);
        else
            orig_prod = NaN;
        end

        ko_indices = find(z0_base == 1 & z_multi < 0.5);
        act_indices = find(z0_base == 0 & z_multi > 0.5);
        all_perturb_indices = [ko_indices(:); act_indices(:)];

        if isempty(all_perturb_indices)
            continue;
        end

        fprintf('Analyzing %s (%d total perturbations)...\n', strat_name, numel(all_perturb_indices));

        for p = 1:numel(all_perturb_indices)
            tf_idx = all_perturb_indices(p);
            tf_name = TF_list{tf_idx};

            z_sKO = z0_base;
            if z_multi(tf_idx) < 0.5
                z_sKO(tf_idx) = 0;
                perturb_type = 'Single_KO';
            else
                z_sKO(tf_idx) = 1;
                perturb_type = 'Single_Activation';
            end

            TF_activity_sKO = TF_exp .* z_sKO;
            y_sKO = beta_matrix * TF_activity_sKO;
            y_sKO = gene_to_rxn_map * y_sKO;
            y_sKO = y_sKO ./ tau;
            y_sKO(isnan(y_sKO)) = 1;
            y_sKO(~isfinite(y_sKO)) = 1;

            A_sKO = sparse(0, nVars_sKO);
            rhs_sKO = [];
            sen_sKO = '';

            A_sKO = [A_sKO; sparse(modelWT_FBA.S), sparse(nMets_sKO, 2*nRxns_sKO)];
            rhs_sKO = [rhs_sKO; zeros(nMets_sKO,1)];
            sen_sKO = [sen_sKO; repmat('=', nMets_sKO, 1)];

            for i = 1:nRxns_sKO
                y_i = y_sKO(i);
                if any(gene_to_rxn_map(i,:)) && y_i < 0.99
                    row_lb = zeros(1, nVars_sKO);
                    row_lb(v0_sKO + i) = 1;
                    row_lb(a0_sKO + i) = 1;
                    A_sKO = [A_sKO; row_lb]; %#ok<AGROW>
                    rhs_sKO = [rhs_sKO; modelWT_FBA.lb(i) * y_i]; %#ok<AGROW>
                    sen_sKO = [sen_sKO; '>']; %#ok<AGROW>

                    row_ub = zeros(1, nVars_sKO);
                    row_ub(v0_sKO + i) = 1;
                    row_ub(b0_sKO + i) = -1;
                    A_sKO = [A_sKO; row_ub]; %#ok<AGROW>
                    rhs_sKO = [rhs_sKO; modelWT_FBA.ub(i) * y_i]; %#ok<AGROW>
                    sen_sKO = [sen_sKO; '<']; %#ok<AGROW>
                end
            end

            lb_sKO = [modelWT_FBA.lb; zeros(2*nRxns_sKO,1)];
            ub_sKO = [modelWT_FBA.ub; inf(2*nRxns_sKO,1)];
            lb_sKO(v0_sKO + biomass_idx) = biomassWT_FBA * flexBiomass;
            ub_sKO(v0_sKO + biomass_idx) = 1000;

            model_sKO = struct();
            model_sKO.A = sparse(A_sKO);
            model_sKO.rhs = rhs_sKO;
            model_sKO.sense = sen_sKO;
            model_sKO.lb = lb_sKO;
            model_sKO.ub = ub_sKO;
            model_sKO.vtype = vtype_sKO;
            model_sKO.modelsense = 'max';

            weight_prod = max(weightMin, kappa_val * weightMultiplier);
            model_sKO.obj = [zeros(nRxns_sKO,1); -kappa_val*ones(nRxns_sKO,1); -kappa_val*ones(nRxns_sKO,1)];
            model_sKO.obj(v0_sKO + prod_idx) = weight_prod;

            gparams = struct();
            gparams.OutputFlag = 0;
            result_sKO = gurobi(model_sKO, gparams);

            if isfield(result_sKO, 'status') && strcmp(result_sKO.status, 'OPTIMAL')
                v_sol_sKO = result_sKO.x(1:nRxns_sKO);
                sKO_prod = v_sol_sKO(prod_idx);
            else
                sKO_prod = NaN;
            end

            if ~isnan(sKO_prod) && ~isnan(orig_prod) && orig_prod > 1e-6
                prod_drop_pct = 100 * (orig_prod - sKO_prod) / orig_prod;
            else
                prod_drop_pct = 0;
            end

            sKO_results = [sKO_results; table({strat_name}, {tf_name}, {perturb_type}, ...
                orig_prod, sKO_prod, prod_drop_pct, ...
                'VariableNames', emptyVarNames)]; %#ok<AGROW>

            if numel(all_perturb_indices) < nKOsLimit
                fprintf('  [%s] %s: Prod went from %.4f -> %.4f (Drop: %5.1f%%%%) \n', ...
                    perturb_type, tf_name, orig_prod, sKO_prod, prod_drop_pct);
            end
        end
    end
end

function value = getParam(params, fieldName, defaultValue)
    if isstruct(params) && isfield(params, fieldName) && ~isempty(params.(fieldName))
        value = params.(fieldName);
    else
        value = defaultValue;
    end
end

function out = toCellArrayLocal(x)
    if isstring(x)
        out = cellstr(x(:));
    elseif iscellstr(x)
        out = x(:);
    elseif iscell(x)
        out = cellfun(@char, x(:), 'UniformOutput', false);
    else
        out = cellstr(string(x(:)));
    end
end
