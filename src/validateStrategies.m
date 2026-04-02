function final_results = validateStrategies(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, maxProduct, ...
        TF_list, WT_kappa, upsilon, resultsTableFilter, params)

    fprintf('\n--- Validation: predicting growth and product with predicted z_sol ---\n');

    if isempty(resultsTableFilter)
        final_results = table();
        return;
    end

    nSols = height(resultsTableFilter);
    KO_names = cell(nSols,1);
    act_names = cell(nSols,1);
    strategy_val = cell(nSols,1);
    iteration_val = NaN(nSols,1);

    growth_val_vbio = NaN(nSols,1);
    prod_val_vbio   = NaN(nSols,1);
    slack_val_vbio  = NaN(nSols,1);

    growth_val_prod = NaN(nSols,1);
    prod_val_prod   = NaN(nSols,1);
    slack_val_prod  = NaN(nSols,1);

    [nMets, nRxns] = size(modelWT_FBA.S);
    nTF = numel(TF_list);

    v0 = 0;
    a0 = v0 + nRxns;
    b0 = a0 + nRxns;
    z0_off = b0 + nRxns;
    nVars = z0_off + nTF;

    A = sparse(0, nVars);
    rhs = [];
    sen = '';

    A = [A; sparse(modelWT_FBA.S), sparse(nMets, nVars - nRxns)];
    rhs = [rhs; zeros(nMets,1)];
    sen = [sen; repmat('=', nMets, 1)];

    for i = 1:nRxns
        if any(upsilon(i,:))
            row = zeros(1, nVars);
            row(v0 + i) = 1;
            row(a0 + i) = 1;
            row(z0_off+1 : z0_off+nTF) = -modelWT_FBA.lb(i) * upsilon(i,:);
            A = [A; row];
            rhs = [rhs; 0];
            sen = [sen; '>'];

            row = zeros(1, nVars);
            row(v0 + i) = 1;
            row(b0 + i) = -1;
            row(z0_off+1 : z0_off+nTF) = modelWT_FBA.ub(i) * upsilon(i,:);
            A = [A; row];
            rhs = [rhs; 0];
            sen = [sen; '<'];
        end
    end

    gparams = struct();
    gparams.OutputFlag = 0;

    for z = 1:nSols
        KO_names{z} = resultsTableFilter.Knockouts{z,1};
        act_names{z} = resultsTableFilter.Activated{z,1};
        strategy_val{z} = resultsTableFilter.Strategy{z,1};
        iteration_val(z) = resultsTableFilter.Iteration(z,1);

        z_val = resultsTableFilter.z_sol{z,1};
        kappa_val = resultsTableFilter.Kappa(z,1);

        fprintf('Validating engineered solution %d of %d\n', z, nSols);

        lb_val = [modelWT_FBA.lb; zeros(2*nRxns,1); z_val];
        ub_val = [modelWT_FBA.ub; inf(2*nRxns,1); z_val];

        lb_val(v0 + biomass_idx) = resultsTableFilter.Biomass(z);
        lb_val(v0 + prod_idx) = resultsTableFilter.Product(z);

        model_val_bio = struct();
        model_val_bio.A = sparse(A);
        model_val_bio.rhs = rhs;
        model_val_bio.sense = sen;
        model_val_bio.lb = lb_val;
        model_val_bio.ub = ub_val;
        model_val_bio.vtype = [repmat('C', 3*nRxns,1); repmat('B', nTF,1)];
        model_val_bio.modelsense = 'max';
        model_val_bio.obj = zeros(nVars,1);
        model_val_bio.obj(v0 + biomass_idx) = params.mu;
        model_val_bio.obj(a0+1 : b0) = -kappa_val;
        model_val_bio.obj(b0+1 : z0_off) = -kappa_val;

        result_val_bio = gurobi(model_val_bio, gparams);

        if isfield(result_val_bio, 'status') && strcmp(result_val_bio.status, 'OPTIMAL')
            v_sol = result_val_bio.x(1:nRxns);
            alpha_vals = result_val_bio.x(nRxns+1:2*nRxns);
            beta_vals = result_val_bio.x(2*nRxns+1:3*nRxns);

            growth_val_vbio(z) = v_sol(biomass_idx);
            prod_val_vbio(z) = v_sol(prod_idx);
            slack_val_vbio(z) = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

            fprintf('  [Bio max] Biomass = %.4f (%.2f%% WT), Product = %.4f (%.2f%% max), slacks = %d\n', ...
                growth_val_vbio(z), 100*growth_val_vbio(z)/biomassWT_FBA, ...
                prod_val_vbio(z), 100*prod_val_vbio(z)/maxProduct, ...
                slack_val_vbio(z));
        end

        model_val_prod = model_val_bio;
        model_val_prod.obj = zeros(nVars,1);
        model_val_prod.obj(v0 + prod_idx) = 1;
        model_val_prod.obj(a0+1 : b0) = -kappa_val;
        model_val_prod.obj(b0+1 : z0_off) = -kappa_val;

        result_val_prod = gurobi(model_val_prod, gparams);

        if isfield(result_val_prod, 'status') && strcmp(result_val_prod.status, 'OPTIMAL')
            v_sol = result_val_prod.x(1:nRxns);
            alpha_vals = result_val_prod.x(nRxns+1:2*nRxns);
            beta_vals = result_val_prod.x(2*nRxns+1:3*nRxns);

            growth_val_prod(z) = v_sol(biomass_idx);
            prod_val_prod(z) = v_sol(prod_idx);
            slack_val_prod(z) = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

            fprintf('  [Prod max] Biomass = %.4f (%.2f%% WT), Product = %.4f (%.2f%% max), slacks = %d\n', ...
                growth_val_prod(z), 100*growth_val_prod(z)/biomassWT_FBA, ...
                prod_val_prod(z), 100*prod_val_prod(z)/maxProduct, ...
                slack_val_prod(z));
        end
    end

    final_results = table('Size', [nSols, 10], ...
        'VariableTypes', {'cell','double','cell','cell','double','double','double','double','double','double'}, ...
        'VariableNames', {'Strategy','Iteration','Knockouts','Activations', ...
                          'Min_growth','Max_growth','Min_product','Max_product', ...
                          'Slack_growth','Slack_product'});

    final_results.Strategy = strategy_val;
    final_results.Iteration = iteration_val;
    final_results.Knockouts = KO_names;
    final_results.Activations = act_names;
    final_results.Max_growth = growth_val_vbio;
    final_results.Min_product = prod_val_vbio;
    final_results.Min_growth = growth_val_prod;
    final_results.Max_product = prod_val_prod;
    final_results.Slack_growth = slack_val_vbio;
    final_results.Slack_product = slack_val_prod;

    KOconcatName = strings(height(final_results),1);
    ACTconcatName = strings(height(final_results),1);

    for i = 1:height(final_results)
        KOconcatName(i) = joinCellString(final_results.Knockouts{i});
        ACTconcatName(i) = joinCellString(final_results.Activations{i});
    end
    final_results.Knockouts = KOconcatName;
    final_results.Activations = ACTconcatName;

    rxnNames = string(modelWT_FBA.rxns);
    affectedRxnNames = strings(height(final_results),1);

    for i = 1:height(final_results)
        tfList = strings(0,1);
        if final_results.Knockouts(i) ~= ""
            tfList = [tfList; string(strtrim(strsplit(final_results.Knockouts(i), ',')))'];
        end
        if final_results.Activations(i) ~= ""
            tfList = [tfList; string(strtrim(strsplit(final_results.Activations(i), ',')))'];
        end
        tfList = unique(tfList(tfList ~= ""));

        if isempty(tfList)
            continue;
        end

        tf_idx = find(ismember(string(TF_list), tfList));
        affected_rxn_idx = find(any(upsilon(:, tf_idx) ~= 0, 2));

        if ~isempty(affected_rxn_idx)
            affectedRxnNames(i) = strjoin(cellstr(rxnNames(affected_rxn_idx)), ', ');
        end
    end

    final_results.AffectedReactions = affectedRxnNames;
end