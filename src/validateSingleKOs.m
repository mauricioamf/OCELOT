function sKO_results = validateSingleKOs(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, maxProduct, ...
        TF_list, z0, WT_kappa, upsilon, final_results, params)

    fprintf('\n--- Validating Single Knockouts ---\n');

    knockout_strings = string(final_results.Knockouts);
    all_genes = strings(0,1);
    for i = 1:numel(knockout_strings)
        if knockout_strings(i) ~= ""
            split_genes = string(strtrim(strsplit(knockout_strings(i), ',')));
            all_genes = [all_genes; split_genes(:)]; %#ok<AGROW>
        end
    end
    all_genes = unique(all_genes(all_genes ~= ""));
    if isempty(all_genes)
        sKO_results = table();
        return;
    end

    knockout_idx = find(ismember(string(TF_list), all_genes));
    KO_names = cellstr(all_genes);
    nKOs = numel(knockout_idx);

    growth_sKO_vbio  = NaN(nKOs,1);
    prod_sKO_vbio    = NaN(nKOs,1);
    growth_sKO_prod  = NaN(nKOs,1);
    prod_sKO_prod    = NaN(nKOs,1);

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

    for k = 1:nKOs
        tf_idx = knockout_idx(k);
        z_val = z0;
        z_val(tf_idx) = 0;

        if nKOs < params.nKOsLimit
            fprintf('Validating KO of TF %s (%d of %d)\n', KO_names{k}, k, nKOs);
        end

        lb_val = [modelWT_FBA.lb; zeros(2*nRxns,1); z_val];
        ub_val = [modelWT_FBA.ub; inf(2*nRxns,1); z_val];

        lb_val(v0 + biomass_idx) = biomassWT_FBA * params.flexBiomass;
        lb_val(v0 + prod_idx) = maxProduct * params.flexProduct;

        model_bio = struct();
        model_bio.A = sparse(A);
        model_bio.rhs = rhs;
        model_bio.sense = sen;
        model_bio.lb = lb_val;
        model_bio.ub = ub_val;
        model_bio.vtype = [repmat('C', 3*nRxns,1); repmat('B', nTF,1)];
        model_bio.modelsense = 'max';
        model_bio.obj = zeros(nVars,1);
        model_bio.obj(v0 + biomass_idx) = params.mu;
        model_bio.obj(a0+1 : b0) = -WT_kappa;
        model_bio.obj(b0+1 : z0_off) = -WT_kappa;

        result_bio = gurobi(model_bio, gparams);
        if isfield(result_bio, 'status') && strcmp(result_bio.status, 'OPTIMAL')
            v_sol = result_bio.x(1:nRxns);
            growth_sKO_vbio(k) = v_sol(biomass_idx);
            prod_sKO_vbio(k) = v_sol(prod_idx);

            if nKOs < params.nKOsLimit
                fprintf('  [Bio max] Biomass = %.4f, Product = %.4f\n', growth_sKO_vbio(k), prod_sKO_vbio(k));
            end
        end

        model_prod = model_bio;
        model_prod.obj = zeros(nVars,1);
        model_prod.obj(v0 + prod_idx) = 1;
        model_prod.obj(a0+1 : b0) = -WT_kappa;
        model_prod.obj(b0+1 : z0_off) = -WT_kappa;

        result_prod = gurobi(model_prod, gparams);
        if isfield(result_prod, 'status') && strcmp(result_prod.status, 'OPTIMAL')
            v_sol = result_prod.x(1:nRxns);
            growth_sKO_prod(k) = v_sol(biomass_idx);
            prod_sKO_prod(k) = v_sol(prod_idx);

            if nKOs < params.nKOsLimit
                fprintf('  [Prod max] Biomass = %.4f, Product = %.4f\n', growth_sKO_prod(k), prod_sKO_prod(k));
            end
        end
    end

    sKO_results = table('Size', [nKOs, 5], ...
        'VariableTypes', {'cell','double','double','double','double'}, ...
        'VariableNames', {'KO_names', 'Max_growth', 'Min_product', 'Min_growth', 'Max_product'});

    sKO_results.KO_names = KO_names;
    sKO_results.Max_growth = growth_sKO_vbio;
    sKO_results.Min_product = prod_sKO_vbio;
    sKO_results.Min_growth = growth_sKO_prod;
    sKO_results.Max_product = prod_sKO_prod;
end