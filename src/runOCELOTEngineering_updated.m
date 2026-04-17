
function results = runOCELOTEngineering(model, expression_data, beta_df, params)
% runOCELOTEngineering
% Self-contained OCELOT engineering workflow with support for heterotrophic,
% mixotrophic, and autotrophic model preparation.
%
% Required inputs
%   model            COBRA model structure
%   expression_data  table with genes in RowNames and samples in VariableNames
%   beta_df          table with target genes in RowNames and TFs in VariableNames
%   params           struct with OCELOT and model-configuration settings
%
% Required params
%   .selectedReference  WT reference condition
%   .biomassRxn         biomass reaction ID
%   .prodRxn            product exchange / demand reaction ID
%
% Growth-mode params
%   .growthType         'heterotrophic' (default), 'mixotrophic', or 'autotrophic'
%
% Heterotrophic / mixotrophic
%   .mediumRxns         exchange reactions opened in the medium
%   .glucoseRxn         glucose exchange reaction ID
%   .glucoseUptake      default -10
%
% Autotrophic (Synechocystis-style)
%   .photonRxn                     default 'EX_photon_e'
%   .glucoseRxn                    default 'EX_glc__D_e'
%   .hco3Rxn                       default 'EX_hco3_e'
%   .co2Rxn                        default 'EX_co2_e'
%   .h2co3TransportRxn             default 'H2CO3_NAt_syn'
%   .HCO3Uptake                    default -3.7
%   .autotrophicOffBiomassRxns     default ["Ec_biomass_SynHetero","Ec_biomass_SynMixo"]
%   .autotrophicDisableRxns        Synechocystis-specific reactions constrained to zero
%   .fnorRxn                       default 'FNOR'
%
% OCELOT params
%   .thresholdSource     'tf_expression' or 'z_vector'
%   .thresholdSelection  'Growth', 'RMSE', 'Discriminability', or 'Composite'
%   .negativeBetaMode    'zero' (default) or 'keep'
%   .upsilonMode         'wt_ratio' or 'global_max'
%   .maxIterations       default 12
%   .maxPerturbations    default 5
%   .biomassFractionStep2 default 0.5
%   .flexBiomass         default 0.5
%   .flexProduct         default 0.5
%   .validateStrategies  default true
%   .validateSingleKOs   default false
%   .validateLeaveOneOut default false
%   .nKOsLimit           default 20
%   .protectedRxns       reactions excluded from regulatory constraints
%   .additionalBoundChanges struct array with fields .rxn .value .bound
%
% Output
%   results  struct containing WT selection, engineering MILP results,
%            validated strategies, and optional single-KO validations.

    params = applyDefaults(params);

    [avg_expr_data, unique_conditions, selectedWTReference] = preprocessExpressionData(expression_data, params.selectedReference);
    [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureBaseModel(model, params);
    reg = buildRegulatoryContext(modelWT_FBA, avg_expr_data, beta_df, selectedWTReference, params);

    wt = runWTThresholdSelection(modelWT_FBA, biomass_idx, reg, protectedRxnIdx, params);
    [best_threshold, best_idx, best_thr_value] = selectBestThreshold( ...
        wt.rmse_table, wt.threshold_names, wt.threshold_values, params.thresholdSelection);

    WT_kappa = wt.WT_results.(best_threshold).kappa;
    z0 = wt.z_struct.(best_threshold);

    prod_idx = find(strcmp(modelWT_FBA.rxns, params.prodRxn), 1);
    if isempty(prod_idx)
        error('Product reaction "%s" was not found in the model.', string(params.prodRxn));
    end

    modelProd_FBA = modelWT_FBA;
    modelProd_FBA.c(:) = 0;
    modelProd_FBA.c(prod_idx) = 1;
    modelProd_FBA.lb(prod_idx) = -1000;
    modelProd_FBA.ub(prod_idx) = 1000;
    modelProd_FBA.lb(biomass_idx) = biomassWT_FBA * params.biomassFractionStep2;
    solProd_FBA = optimizeCbModel(modelProd_FBA);

    y_WT_base = reg.gene_to_rxn_map * (reg.beta_matrix * (reg.TF_exp .* z0));
    y_WT_base = replaceInfWithFinite(y_WT_base);

    upsilon = buildUpsilon(reg.gene_to_rxn_map, reg.beta_prime, y_WT_base, protectedRxnIdx, params.upsilonMode);

    [resultsTable, rawResults] = runEngineeringMILP(modelWT_FBA, biomass_idx, prod_idx, ...
        biomassWT_FBA, solProd_FBA.f, reg.TF_list, z0, upsilon, WT_kappa, params);

    resultsTableFilter = filterStrategies(resultsTable);

    if params.validateStrategies && ~isempty(resultsTableFilter)
        final_results = validateStrategies(modelWT_FBA, biomass_idx, prod_idx, ...
            biomassWT_FBA, solProd_FBA.f, reg.TF_list, WT_kappa, upsilon, ...
            resultsTableFilter, params);
    else
        final_results = table();
    end

    tau = max(abs(y_WT_base));
    if tau == 0
        tau = 1;
    end

    if params.validateSingleKOs && ~isempty(resultsTableFilter) && ~isempty(final_results)
        sKO_results = validateSingleKOs(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, ...
            reg.TF_list, z0, reg.TF_exp, reg.beta_matrix, reg.gene_to_rxn_map, tau, ...
            resultsTableFilter, final_results, params);
    else
        sKO_results = table();
    end

    if params.validateLeaveOneOut && ~isempty(resultsTableFilter) && ~isempty(final_results)
        loo_results = validateLeaveOneOut(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, ...
            reg.TF_list, z0, reg.TF_exp, reg.beta_matrix, reg.gene_to_rxn_map, tau, ...
            resultsTableFilter, final_results, params);
    else
        loo_results = table();
    end

    results = struct();
    results.paramsUsed            = params;
    results.uniqueConditions      = unique_conditions;
    results.selectedWTReference   = selectedWTReference;
    results.avgExprData           = avg_expr_data;
    results.modelWT_FBA           = modelWT_FBA;
    results.FBAsol                = FBAsol;
    results.biomassWT_FBA         = biomassWT_FBA;
    results.protectedRxns         = protectedRxns;
    results.protectedRxnIdx       = protectedRxnIdx;
    results.solProd_FBA           = solProd_FBA;
    results.WT_results            = wt.WT_results;
    results.rmse_table            = wt.rmse_table;
    results.best_threshold        = best_threshold;
    results.best_threshold_value  = best_thr_value;
    results.best_threshold_index  = best_idx;
    results.WT_kappa              = WT_kappa;
    results.z0                    = z0;
    results.rawResultsCell        = rawResults;
    results.resultsTable          = resultsTable;
    results.filteredResultsTable  = resultsTableFilter;
    results.validatedResultsTable = final_results;
    results.singleKOResults       = sKO_results;
    results.leaveOneOutResults    = loo_results;
    results.regulatoryContext     = reg;
    results.upsilon               = upsilon;
    results.y_WT_base             = y_WT_base;
end

% ========================= LOCAL FUNCTIONS =========================

function params = applyDefaults(params)
    if nargin < 1 || isempty(params), params = struct(); end

    if ~isfield(params, 'growthType') || isempty(params.growthType)
        params.growthType = 'heterotrophic';
    end

    growthType = lower(string(params.growthType));

    defaults.mu                 = 1;
    defaults.maxInfeas          = 3;
    defaults.kappaStart         = 1000;
    defaults.kappaDecay         = 10;
    defaults.removeLastExchange = true;
    defaults.negativeBetaMode   = 'zero';
    defaults.protectedRxns      = strings(0,1);
    defaults.additionalBoundChanges = struct([]);
    defaults.glucoseUptake      = -10;
    defaults.mediumRxns         = strings(0,1);
    defaults.glucoseRxn         = 'EX_glc__D_e';
    defaults.photonRxn          = 'EX_photon_e';
    defaults.hco3Rxn            = 'EX_hco3_e';
    defaults.co2Rxn             = 'EX_co2_e';
    defaults.h2co3TransportRxn  = 'H2CO3_NAt_syn';
    defaults.HCO3Uptake         = -3.7;
    defaults.fnorRxn            = 'FNOR';
    defaults.autotrophicOffBiomassRxns = ["Ec_biomass_SynHetero","Ec_biomass_SynMixo"];
    defaults.autotrophicDisableRxns = [ ...
        "CBFC2ub","CBFC2pb","CYO1b_syn","CYO1bpp_syn","PSI_2a", ...
        "NDH2_1p","NDH2_syn","NDH1_1p","CYO1b2pp_syn","CBFCpb", ...
        "CYTBDpp","SUCDyy_syn","FDH6pp","G3PDap","PROD5p", ...
        "FDH6","PROD5u","G3PDau"];
    defaults.thresholdSource = '';
    defaults.thresholdSelection = '';
    defaults.kappaMin = [];
    defaults.upsilonMode = '';
    defaults.maxIterations = 12;
    defaults.maxPerturbations = 5;
    defaults.biomassFractionStep2 = 0.5;
    defaults.flexBiomass = 0.5;
    defaults.flexProduct = 0.5;
    defaults.timeLimit = 600;
    defaults.mipGap = 0.01;
    defaults.mipFocus = 1;
    defaults.validateStrategies = true;
    defaults.validateSingleKOs = false;
    defaults.nKOsLimit = 20;

    fn = fieldnames(defaults);
    for i = 1:numel(fn)
        if ~isfield(params, fn{i}) || isempty(params.(fn{i}))
            params.(fn{i}) = defaults.(fn{i});
        end
    end

    if strlength(string(params.thresholdSource)) == 0
        if growthType == "autotrophic"
            params.thresholdSource = 'z_vector';
        else
            params.thresholdSource = 'tf_expression';
        end
    end

    if strlength(string(params.thresholdSelection)) == 0
        if growthType == "autotrophic"
            params.thresholdSelection = 'RMSE';
        else
            params.thresholdSelection = 'Composite';
        end
    end

    if isempty(params.kappaMin)
        if growthType == "autotrophic"
            params.kappaMin = 1e-12;
        else
            params.kappaMin = params.kappaStart;
        end
    end

    if strlength(string(params.upsilonMode)) == 0
        if growthType == "autotrophic"
            params.upsilonMode = 'global_max';
        else
            params.upsilonMode = 'wt_ratio';
        end
    end

    required = {'selectedReference','biomassRxn','prodRxn'};
    for i = 1:numel(required)
        if ~isfield(params, required{i}) || isempty(params.(required{i}))
            error('params.%s is required.', required{i});
        end
    end

    if growthType ~= "autotrophic"
        if ~isfield(params, 'mediumRxns') || isempty(params.mediumRxns)
            error('params.mediumRxns is required for %s growth.', char(growthType));
        end
    end
end

function [avg_expr_data, unique_conditions, selectedWTReference] = preprocessExpressionData(expression_data, selectedReference)
    all_samples = expression_data.Properties.VariableNames;
    condition_names = regexprep(all_samples, '__\d+$', '');

    for i = 1:numel(condition_names)
        if ~isempty(condition_names{i}) && isstrprop(condition_names{i}(1), 'digit')
            condition_names{i} = ['X' condition_names{i}];
        end
    end

    unique_conditions = unique(condition_names, 'stable');
    wt_idx = find(contains(unique_conditions, selectedReference), 1, 'first');
    if isempty(wt_idx)
        error('WT reference "%s" not found in expression_data.', string(selectedReference));
    end
    selectedWTReference = unique_conditions{wt_idx};

    avg_expr_data = array2table( ...
        zeros(size(expression_data,1), numel(unique_conditions)), ...
        'RowNames', expression_data.Properties.RowNames, ...
        'VariableNames', unique_conditions);

    for c = 1:numel(unique_conditions)
        cond = unique_conditions{c};
        idx = strcmp(condition_names, cond);
        avg_expr_data{:,c} = mean(expression_data{:,idx}, 2, 'omitnan');
    end
end

function [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureBaseModel(model, params)
    growthType = lower(string(params.growthType));

    switch growthType
        case "autotrophic"
            [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureAutotrophicModel(model, params);
        otherwise
            modelWT_FBA = model;

            exchangeIndex = find(contains(lower(string(modelWT_FBA.rxnNames)), "exchange"));
            exchangeIDs = string(modelWT_FBA.rxns(exchangeIndex));
            if params.removeLastExchange && ~isempty(exchangeIDs)
                exchangeIDs(end) = [];
            end

            biomass_idx = find(strcmp(modelWT_FBA.rxns, params.biomassRxn), 1);
            if isempty(biomass_idx)
                error('Biomass reaction "%s" not found in the model.', string(params.biomassRxn));
            end

            if ~isempty(exchangeIDs)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, cellstr(exchangeIDs), 0, 'l');
            end
            if ~isempty(params.mediumRxns)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, cellstr(toStringArray(params.mediumRxns)), -1000, 'l');
            end
            modelWT_FBA = changeObjective(modelWT_FBA, modelWT_FBA.rxns{biomass_idx});

            if ~isempty(params.glucoseRxn)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, char(params.glucoseRxn), params.glucoseUptake, 'b');
            end

            modelWT_FBA = applyBoundChangeList(modelWT_FBA, params.additionalBoundChanges);

            FBAsol = optimizeCbModel(modelWT_FBA);
            biomassWT_FBA = FBAsol.f;

            protectedRxns = unique([toStringArray(params.protectedRxns); strings(0,1)]);
            protectedRxnIdx = find(ismember(string(modelWT_FBA.rxns), protectedRxns));
    end
end

function [modelAuto, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureAutotrophicModel(model, params)
    modelAuto = model;

    biomass_idx = find(strcmp(modelAuto.rxns, params.biomassRxn), 1);
    if isempty(biomass_idx)
        error('Biomass reaction "%s" not found in the model.', string(params.biomassRxn));
    end

    modelAuto = changeObjective(modelAuto, char(params.biomassRxn));

    offBiomass = toStringArray(params.autotrophicOffBiomassRxns);
    for i = 1:numel(offBiomass)
        if any(strcmp(modelAuto.rxns, offBiomass(i)))
            modelAuto = changeRxnBounds(modelAuto, char(offBiomass(i)), 0, 'b');
        end
    end

    if ~isempty(params.photonRxn) && any(strcmp(modelAuto.rxns, params.photonRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.photonRxn), -1000, 'l');
    end
    if ~isempty(params.glucoseRxn) && any(strcmp(modelAuto.rxns, params.glucoseRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.glucoseRxn), 0, 'b');
    end
    if ~isempty(params.hco3Rxn) && any(strcmp(modelAuto.rxns, params.hco3Rxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.hco3Rxn), params.HCO3Uptake, 'l');
    end
    if ~isempty(params.co2Rxn) && any(strcmp(modelAuto.rxns, params.co2Rxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.co2Rxn), 0, 'b');
    end
    if ~isempty(params.h2co3TransportRxn) && any(strcmp(modelAuto.rxns, params.h2co3TransportRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.h2co3TransportRxn), 0, 'u');
    end

    disableRxns = toStringArray(params.autotrophicDisableRxns);
    for i = 1:numel(disableRxns)
        if any(strcmp(modelAuto.rxns, disableRxns(i)))
            modelAuto = changeRxnBounds(modelAuto, char(disableRxns(i)), 0, 'b');
        end
    end

    if ~isempty(params.fnorRxn) && any(strcmp(modelAuto.rxns, params.fnorRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.fnorRxn), 0, 'l');
    end

    modelAuto = applyBoundChangeList(modelAuto, params.additionalBoundChanges);

    solAuto1 = optimizeCbModel(modelAuto, 'max');
    if isempty(solAuto1) || ~isfield(solAuto1, 'f') || isempty(solAuto1.f)
        error('Autotrophic WT optimization failed at the first biomass step.');
    end

    modelAuto2 = changeObjective(modelAuto, char(params.photonRxn), 1);
    modelAuto2 = changeRxnBounds(modelAuto2, char(params.biomassRxn), solAuto1.f, 'b');
    solAuto2 = optimizeCbModel(modelAuto2, 'max');

    if isempty(solAuto2) || ~isfield(solAuto2, 'f') || isempty(solAuto2.f)
        error('Autotrophic photon minimization failed.');
    end

    modelAuto = changeRxnBounds(modelAuto, char(params.photonRxn), solAuto2.f, 'l');
    solAuto3 = optimizeCbModel(modelAuto, 'max');

    if isempty(solAuto3) || ~isfield(solAuto3, 'x') || isempty(solAuto3.x)
        error('Autotrophic WT optimization failed after fixing photon uptake.');
    end

    biomassWT_FBA = solAuto3.x(biomass_idx);
    FBAsol = solAuto3;

    protectedRxns = unique([ ...
        toStringArray(params.protectedRxns); ...
        toStringArray(params.autotrophicOffBiomassRxns); ...
        string(params.photonRxn); string(params.glucoseRxn); string(params.hco3Rxn); ...
        string(params.co2Rxn); string(params.h2co3TransportRxn); string(params.fnorRxn); ...
        toStringArray(params.autotrophicDisableRxns)]);

    protectedRxns = protectedRxns(strlength(protectedRxns) > 0);
    protectedRxnIdx = find(ismember(string(modelAuto.rxns), protectedRxns));
end

function model = applyBoundChangeList(model, changes)
    if isempty(changes)
        return;
    end

    if ~isstruct(changes)
        error('additionalBoundChanges must be a struct array with fields rxn, value, and bound.');
    end

    for i = 1:numel(changes)
        if ~isfield(changes(i), 'rxn') || ~isfield(changes(i), 'value') || ~isfield(changes(i), 'bound')
            error('Each bound-change entry must contain fields rxn, value, and bound.');
        end
        if any(strcmp(model.rxns, changes(i).rxn))
            model = changeRxnBounds(model, char(changes(i).rxn), changes(i).value, changes(i).bound);
        end
    end
end

function reg = buildRegulatoryContext(modelWT_FBA, avg_expr_data, beta_df, selectedWTReference, params)
    TF_list = cellstr(string(beta_df.Properties.VariableNames(:)));
    target_genes = beta_df.Properties.RowNames;

    missingTFs = setdiff(TF_list, avg_expr_data.Properties.RowNames);
    if ~isempty(missingTFs)
        error('The following TFs are missing from expression_data.RowNames: %s', strjoin(missingTFs, ', '));
    end

    missingTargets = setdiff(target_genes, avg_expr_data.Properties.RowNames);
    if ~isempty(missingTargets)
        error('The following target genes are missing from expression_data.RowNames: %s', strjoin(missingTargets, ', '));
    end

    TF_exp = avg_expr_data{TF_list, selectedWTReference};
    target_exp = avg_expr_data{target_genes, selectedWTReference};

    [~, targetGeneIdx_inModel] = ismember(target_genes, modelWT_FBA.genes);
    validTargets = targetGeneIdx_inModel > 0;

    targetGeneIdx_inModel = targetGeneIdx_inModel(validTargets);
    beta_df_valid = beta_df(validTargets, :);
    target_exp_valid = target_exp(validTargets);

    modelWT_FBA = buildRxnGeneMat(modelWT_FBA);
    gene_to_rxn_map = modelWT_FBA.rxnGeneMat(:, targetGeneIdx_inModel);

    beta_matrix = beta_df_valid{:,:};
    if strcmpi(params.negativeBetaMode, 'zero')
        beta_matrix(beta_matrix < 0) = 0;
    end
    beta_prime = beta_matrix .* TF_exp';
    z_vector = pinv(beta_prime) * target_exp_valid;

    threshold_values = chooseThresholdValues(TF_exp, z_vector, params.thresholdSource);

    reg = struct();
    reg.TF_list          = TF_list;
    reg.TF_exp           = TF_exp;
    reg.target_exp       = target_exp;
    reg.target_exp_valid = target_exp_valid;
    reg.target_genes     = target_genes;
    reg.validTargets     = validTargets;
    reg.beta_df_valid    = beta_df_valid;
    reg.beta_matrix      = beta_matrix;
    reg.beta_prime       = beta_prime;
    reg.gene_to_rxn_map  = gene_to_rxn_map;
    reg.z_vector         = z_vector;
    reg.threshold_values = threshold_values;
end

function threshold_values = chooseThresholdValues(TF_exp, z_vector, thresholdSource)
    switch lower(string(thresholdSource))
        case "z_vector"
            threshold_values = unique(sort(z_vector(~isnan(z_vector) & isfinite(z_vector))));
            if isempty(threshold_values), threshold_values = 0; end
        otherwise
            threshold_values = unique(sort(TF_exp(~isnan(TF_exp) & isfinite(TF_exp))));
            if isempty(threshold_values), threshold_values = 0; end
    end
end

function state = buildThresholdState(reg, thresholdValue, thresholdSource)
    switch lower(string(thresholdSource))
        case "z_vector"
            z_bin = double(reg.z_vector >= thresholdValue);
            z_cont = reg.z_vector;
            z_cont(z_cont <= thresholdValue) = 0;
            TF_activity = reg.TF_exp .* z_bin;
            y_tgt = reg.beta_matrix * TF_activity;
            target_exp_hat = reg.beta_prime * z_cont;
        otherwise
            z_bin = double(reg.TF_exp >= thresholdValue);
            z_cont = reg.TF_exp;
            z_cont(z_cont <= thresholdValue) = 0;
            TF_activity = reg.TF_exp .* z_bin;
            y_tgt = reg.beta_matrix * TF_activity;
            target_exp_hat = y_tgt;
    end

    state = struct();
    state.z_bin = z_bin;
    state.z_cont = z_cont;
    state.y_tgt = y_tgt;
    state.target_exp_hat = target_exp_hat;
end

function wt = runWTThresholdSelection(modelWT_FBA, biomass_idx, reg, protectedRxnIdx, params)
    y_all_struct = struct();
    z_struct = struct();
    target_expr_struct = struct();
    threshold_names = cell(numel(reg.threshold_values), 1);

    for t = 1:numel(reg.threshold_values)
        thr = reg.threshold_values(t);
        state = buildThresholdState(reg, thr, params.thresholdSource);

        y_rxn = reg.gene_to_rxn_map * state.y_tgt;
        y_rxn = replaceInfWithFinite(y_rxn);

        T = array2table(y_rxn, 'RowNames', modelWT_FBA.rxns, 'VariableNames', {'WT'});
        tau = max(abs(T{:,'WT'}));
        if tau == 0, tau = 1; end
        T = array2table(T{:,:} ./ tau, 'RowNames', modelWT_FBA.rxns, 'VariableNames', {'WT'});

        name = matlab.lang.makeValidName(sprintf('thr_expr_%g', thr));
        threshold_names{t} = name;
        y_all_struct.(name) = T;
        z_struct.(name) = state.z_bin;
        target_expr_struct.(name) = state.target_exp_hat;
    end

    WT_results = struct();
    infeas_count = 0;

    for tr = 1:numel(threshold_names)
        threshold_name = threshold_names{tr};
        fprintf('\n--- WT PROM for threshold: %s ---\n', threshold_name);

        y_all_loop = y_all_struct.(threshold_name);
        yWT = y_all_loop{:,'WT'};
        target_exp_hat = target_expr_struct.(threshold_name);

        try
            [biomass, gp_WT, result_WT, kappa] = solvePromWithScaledBounds( ...
                modelWT_FBA, biomass_idx, yWT, reg.gene_to_rxn_map, protectedRxnIdx, ...
                params.mu, params.kappaStart, params.kappaMin, params.kappaDecay);

            R = struct();
            R.biomass = biomass;
            R.model = gp_WT;
            R.result = result_WT;
            R.kappa = kappa;
            R.target_expr = array2table(target_exp_hat, ...
                'VariableNames', {threshold_name}, ...
                'RowNames', reg.beta_df_valid.Properties.RowNames);
            WT_results.(threshold_name) = R;
            infeas_count = 0;
        catch
            infeas_count = infeas_count + 1;
            warning('No feasible solution found for %s.', threshold_name);

            R = struct();
            R.biomass = NaN;
            R.model = [];
            R.result = [];
            R.kappa = NaN;
            R.target_expr = array2table(target_exp_hat, ...
                'VariableNames', {threshold_name}, ...
                'RowNames', reg.beta_df_valid.Properties.RowNames);
            WT_results.(threshold_name) = R;

            if infeas_count >= params.maxInfeas
                fprintf('\nStopping WT simulations after %d consecutive infeasible thresholds.\n', params.maxInfeas);
                break;
            end
        end
    end

    simulatedThresholds = fieldnames(WT_results);
    rmse_table = table('Size', [numel(simulatedThresholds), 5], ...
        'VariableTypes', {'string','double','double','double','double'}, ...
        'VariableNames', {'Threshold','RMSE','R2','Growth','Discriminability'});

    for th = 1:numel(simulatedThresholds)
        threshold_name = simulatedThresholds{th};
        predicted_tbl = WT_results.(threshold_name).target_expr;
        predicted_vec = predicted_tbl{:,1};
        target_exp_valid = reg.target_exp_valid;

        predicted_vec = log(abs(predicted_vec) + eps);
        target_exp_valid = log(abs(target_exp_valid) + eps);

        residual = target_exp_valid - predicted_vec;
        goodIdx = isfinite(residual);
        predicted_vec = predicted_vec(goodIdx);
        target_exp_valid = target_exp_valid(goodIdx);

        residual = target_exp_valid - predicted_vec;
        rmse_value = sqrt(mean(residual.^2));
        ss_res = sum((predicted_vec - target_exp_valid).^2);
        ss_tot = sum((target_exp_valid - mean(target_exp_valid)).^2);
        r2 = 1 - (ss_res / (ss_tot + eps));

        state = buildThresholdState(reg, reg.threshold_values(th), params.thresholdSource);
        y_base = reg.gene_to_rxn_map * state.y_tgt;
        delta_norms = zeros(numel(reg.TF_list), 1);
        for i = 1:numel(reg.TF_list)
            if state.z_bin(i) == 0
                continue;
            end
            zko = state.z_bin;
            zko(i) = 0;
            y_ko = reg.gene_to_rxn_map * (reg.beta_matrix * (reg.TF_exp .* zko));
            delta_norms(i) = norm(y_base - y_ko);
        end
        active_deltas = delta_norms(state.z_bin == 1);
        if isempty(active_deltas) || mean(active_deltas) <= 0
            discriminability = 0;
        else
            discriminability = std(active_deltas) / mean(active_deltas);
        end

        rmse_table.Threshold(th) = string(threshold_name);
        rmse_table.RMSE(th) = rmse_value;
        rmse_table.R2(th) = r2;
        rmse_table.Growth(th) = WT_results.(threshold_name).biomass;
        rmse_table.Discriminability(th) = discriminability;
    end

    wt = struct();
    wt.WT_results = WT_results;
    wt.rmse_table = rmse_table;
    wt.threshold_names = simulatedThresholds;
    wt.z_struct = z_struct;
    wt.threshold_values = reg.threshold_values(1:numel(threshold_names));
end

function [best_threshold, best_idx, best_thr_value] = selectBestThreshold(rmse_table, threshold_names, threshold_values, thresholdSelection)
    valid_idx = ~isnan(rmse_table.RMSE) & ~isnan(rmse_table.Growth) & rmse_table.Growth > 0;
    if ~any(valid_idx)
        error('No valid threshold could be selected.');
    end

    switch lower(string(thresholdSelection))
        case "growth"
            [~, best_local_idx] = max(rmse_table.Growth(valid_idx));
        case "rmse"
            [~, best_local_idx] = min(rmse_table.RMSE(valid_idx));
        case "discriminability"
            [~, best_local_idx] = max(rmse_table.Discriminability(valid_idx));
        otherwise
            rmse_valid   = rmse_table.RMSE(valid_idx);
            disc_valid   = rmse_table.Discriminability(valid_idx);
            growth_valid = rmse_table.Growth(valid_idx);
            r2_valid     = rmse_table.R2(valid_idx);

            rmse_norm   = 1 - (rmse_valid - min(rmse_valid)) ./ (max(rmse_valid) - min(rmse_valid) + eps);
            disc_norm   =     (disc_valid - min(disc_valid)) ./ (max(disc_valid) - min(disc_valid) + eps);
            growth_norm =     (growth_valid - min(growth_valid)) ./ (max(growth_valid) - min(growth_valid) + eps);
            r2_norm     = max(r2_valid, 0);

            composite = 0.4 * disc_norm + 0.3 * r2_norm + 0.2 * rmse_norm + 0.1 * growth_norm;
            [~, best_local_idx] = max(composite);
    end

    best_idx_all = find(valid_idx);
    best_idx = best_idx_all(best_local_idx);
    best_threshold = threshold_names{best_idx};
    best_thr_value = threshold_values(best_idx);

    fprintf('\nSelected threshold: %s\n', best_threshold);
    fprintf('  RMSE=%.4f | R2=%.4f | Growth=%.4f | Discriminability=%.4f\n', ...
        rmse_table.RMSE(best_idx), ...
        rmse_table.R2(best_idx), ...
        rmse_table.Growth(best_idx), ...
        rmse_table.Discriminability(best_idx));
end

function upsilon = buildUpsilon(gene_to_rxn_map, beta_prime, y_WT_base, protectedRxnIdx, upsilonMode)
    upsilon = gene_to_rxn_map * beta_prime;

    switch lower(string(upsilonMode))
        case "global_max"
            denom = max(abs(upsilon(:)));
            if denom > 0
                upsilon = upsilon ./ denom;
            end
        otherwise
            for i = 1:size(upsilon,1)
                if y_WT_base(i) > 1e-9
                    upsilon(i,:) = upsilon(i,:) ./ y_WT_base(i);
                else
                    upsilon(i,:) = 0;
                end
            end
    end

    upsilon = max(upsilon, 0);
    upsilon = min(upsilon, 1);

    if ~isempty(protectedRxnIdx)
        upsilon(protectedRxnIdx,:) = 0;
    end
end

function [resultsTable, results_all] = runEngineeringMILP(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, maxProduct, TF_list, z0, upsilon, kappa, params)
    [nMets, nRxns] = size(modelWT_FBA.S);
    nTF = numel(TF_list);

    v0     = 0;
    a0     = v0 + nRxns;
    b0     = a0 + nRxns;
    z0_off = b0 + nRxns;
    nVars  = z0_off + nTF;

    vtype = [repmat('C', nRxns,1); repmat('C', nRxns,1); repmat('C', nRxns,1); repmat('B', nTF,1)];

    gparams = struct();
    gparams.OutputFlag = 0;
    gparams.TimeLimit = params.timeLimit;
    gparams.MIPGap = params.mipGap;
    gparams.MIPFocus = params.mipFocus;

    flexBiomass = params.flexBiomass(:)';
    flexProduct = params.flexProduct(:)';

    if numel(flexBiomass) == 1 && numel(flexProduct) > 1
        flexBiomass = repmat(flexBiomass, size(flexProduct));
    elseif numel(flexProduct) == 1 && numel(flexBiomass) > 1
        flexProduct = repmat(flexProduct, size(flexBiomass));
    elseif numel(flexBiomass) ~= numel(flexProduct)
        error('flexBiomass and flexProduct must have the same length, or one must be scalar.');
    end

    results_all = {};
    integer_cuts = {};

    for idxFlex = 1:numel(flexBiomass)
        fb = flexBiomass(idxFlex);
        fp = flexProduct(idxFlex);

        min_biomass_step3 = biomassWT_FBA * fb;
        min_product_step3 = maxProduct * fp;

        lb_base = modelWT_FBA.lb;
        ub_base = modelWT_FBA.ub;
        lb_base(biomass_idx) = min_biomass_step3;
        lb_base(prod_idx) = min_product_step3;

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
                row(z0_off+1 : z0_off+nTF) = -lb_base(i) * upsilon(i,:);
                A = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '>'];

                row = zeros(1, nVars);
                row(v0 + i) = 1;
                row(b0 + i) = -1;
                row(z0_off+1 : z0_off+nTF) = ub_base(i) * upsilon(i,:);
                A = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '<'];
            end
        end

        perturb_row = zeros(1, nVars);
        for j = 1:nTF
            if z0(j) == 1
                perturb_row(z0_off + j) = -1;
            else
                perturb_row(z0_off + j) = 1;
            end
        end
        rhs_perturb = params.maxPerturbations - sum(z0 == 1);

        A = [A; perturb_row];
        rhs = [rhs; rhs_perturb];
        sen = [sen; '<'];

        lb = [lb_base; zeros(2*nRxns,1); zeros(nTF,1)];
        ub = [ub_base; inf(2*nRxns,1); ones(nTF,1)];

        for iteration = 1:params.maxIterations
            A_with_cuts = A;
            rhs_with_cuts = rhs;
            sen_with_cuts = sen;

            for cut_idx = 1:numel(integer_cuts)
                cut = integer_cuts{cut_idx};
                row = zeros(1, nVars);
                for j = 1:nTF
                    if cut.z_pattern(j) == 1
                        row(z0_off + j) = 1;
                    else
                        row(z0_off + j) = -1;
                    end
                end
                A_with_cuts = [A_with_cuts; row];
                rhs_with_cuts = [rhs_with_cuts; cut.sum_value - 2];
                sen_with_cuts = [sen_with_cuts; '<'];
            end

            obj = zeros(nVars,1);
            for j = 1:nTF
                if z0(j) == 1
                    obj(z0_off + j) = 1;
                else
                    obj(z0_off + j) = -1;
                end
            end
            obj(a0+1 : b0) = -kappa;
            obj(b0+1 : z0_off) = -kappa;

            model_step3 = struct();
            model_step3.A = sparse(A_with_cuts);
            model_step3.rhs = rhs_with_cuts;
            model_step3.sense = sen_with_cuts;
            model_step3.lb = lb;
            model_step3.ub = ub;
            model_step3.vtype = vtype;
            model_step3.obj = obj;
            model_step3.modelsense = 'max';

            result_step3 = gurobi(model_step3, gparams);

            if isfield(result_step3, 'status') && strcmp(result_step3.status, 'OPTIMAL')
                z_sol = result_step3.x(z0_off+1 : z0_off+nTF);
                v_sol = result_step3.x(1:nRxns);
                alpha_vals = result_step3.x(a0+1 : b0);
                beta_vals = result_step3.x(b0+1 : z0_off);

                knockout_idx = find(z0 == 1 & z_sol < 0.5);
                activated_idx = find(z0 == 0 & z_sol > 0.5);

                knocked_out_TFs = TF_list(knockout_idx);
                activated_TFs = TF_list(activated_idx);

                biomass_eng = v_sol(biomass_idx);
                product_eng = v_sol(prod_idx);
                num_slacks_used = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

                results_all(end+1,:) = { ...
                    fb, fp, biomass_eng, product_eng, ...
                    knocked_out_TFs, activated_TFs, ...
                    numel(knocked_out_TFs), numel(activated_TFs), ...
                    num_slacks_used, kappa, z0, z_sol, iteration, 'Original'}; %#ok<AGROW>

                cut_pattern = struct();
                cut_pattern.z_pattern = round(z_sol);
                cut_pattern.sum_value = sum(cut_pattern.z_pattern);

                is_new_pattern = true;
                for cut_idx = 1:numel(integer_cuts)
                    if isequal(integer_cuts{cut_idx}.z_pattern, cut_pattern.z_pattern)
                        is_new_pattern = false;
                        break;
                    end
                end
                if is_new_pattern
                    integer_cuts{end+1} = cut_pattern; %#ok<AGROW>
                end
            else
                if iteration == 1
                    results_all(end+1,:) = {fb, fp, NaN, NaN, {}, {}, 0, 0, NaN, kappa, z0, NaN, 0, 'Original'}; %#ok<AGROW>
                end
                break;
            end
        end
    end

    if isempty(results_all)
        resultsTable = table();
        return;
    end

    resultsTable = cell2table(results_all, ...
        'VariableNames', {'BiomassFlex','ProductFlex','Biomass','Product', ...
                          'Knockouts','Activated','NumKO','NumActivated', ...
                          'SlacksUsed','Kappa','z0','z_sol','Iteration','Strategy'});

    resultsTable = resultsTable(:, {'Strategy','Iteration','BiomassFlex','ProductFlex','Kappa','Biomass','Product', ...
                                    'Knockouts','Activated','NumKO','NumActivated','SlacksUsed','z0','z_sol'});
end

function resultsTableFilter = filterStrategies(resultsTable)
    if isempty(resultsTable)
        resultsTableFilter = table();
        return;
    end

    numericVars = vartype('numeric');
    numericData = table2array(resultsTable(:, numericVars));
    rowHasNaN = any(isnan(numericData), 2);
    rowNoPerturb = (resultsTable.NumKO == 0) & (resultsTable.NumActivated == 0);
    rowAllKOorAllAct = (resultsTable.NumKO == numel(resultsTable.z0{1})) | (resultsTable.NumActivated == numel(resultsTable.z0{1}));
    rowZallZeroes = cellfun(@(v) isempty(v) || all(v == 0), resultsTable.z_sol);

    badRows = rowHasNaN | rowNoPerturb | rowAllKOorAllAct | rowZallZeroes;
    fprintf('Filtering %d / %d rows (NaN / no-perturb / all-TF).\n', sum(badRows), height(resultsTable));
    resultsTableFilter = resultsTable(~badRows, :);
end

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

function out = joinCellString(c)
    if isempty(c)
        out = "";
    else
        out = string(strjoin(cellstr(string(c)), ', '));
    end
end

function [biomass, gp, result, kappa] = solvePromWithScaledBounds(model, biomass_idx, yvec, gene_to_rxn_map, protectedRxnIdx, mu, kappaStart, kappaMin, kappaDecay)
    S = sparse(model.S);
    [m, n] = size(S);
    lb_base = model.lb;
    ub_base = model.ub;

    lbp = lb_base;
    ubp = ub_base;

    isProtected = false(n,1);
    isProtected(protectedRxnIdx) = true;

    yvec(~isfinite(yvec)) = 1;

    for i = 1:n
        if isProtected(i)
            continue;
        end
        if ~any(gene_to_rxn_map(i,:))
            continue;
        end

        y = yvec(i);

        if abs(lb_base(i)) > 1e-6 && abs(ub_base(i)) > 1e-6
            lbp(i) = max(lb_base(i) * y, -1000);
            ubp(i) = min(ub_base(i) * y,  1000);
        elseif lb_base(i) < -1e-6
            lbp(i) = max(lb_base(i) * y, -1000);
        elseif ub_base(i) > 1e-6
            ubp(i) = min(ub_base(i) * y,  1000);
        end
    end

    lbp(~isfinite(lbp)) = lb_base(~isfinite(lbp));
    ubp(~isfinite(ubp)) = ub_base(~isfinite(ubp));

    infeasible_idx = find(lbp > ubp);
    for k = infeasible_idx'
        lbp(k) = min(lbp(k), ubp(k) - 1e-6);
        ubp(k) = max(lbp(k), ubp(k) + 1e-6);
    end

    found_solution = false;
    gp = [];
    result = [];
    biomass = NaN;
    kappa = NaN;

    kappaLoop = kappaStart;
    while kappaLoop >= kappaMin
        gp.A = sparse([ ...
            S,          sparse(m,n),      sparse(m,n); ...
           -speye(n),  -speye(n),         sparse(n,n); ...
            speye(n),   sparse(n,n),     -speye(n)]);

        gp.rhs        = [zeros(m,1); -lbp; ubp];
        gp.sense      = [repmat('=',m,1); repmat('<',n,1); repmat('<',n,1)];
        gp.modelsense = 'min';
        gp.obj        = [zeros(n,1); kappaLoop*ones(n,1); kappaLoop*ones(n,1)];
        gp.obj(biomass_idx) = -mu;
        gp.lb         = [lb_base; zeros(2*n,1)];
        gp.ub         = [ub_base; inf(2*n,1)];

        gparams = struct();
        gparams.OutputFlag = 0;
        resultLoop = gurobi(gp, gparams);

        if isfield(resultLoop, 'x')
            v_sol = resultLoop.x(1:n);
            biomassLoop = v_sol(biomass_idx);

            if biomassLoop > 1e-6
                fprintf('→ predicted biomass = %.4f, kappa = %.3e\n', biomassLoop, kappaLoop);
                biomass = biomassLoop;
                result = resultLoop;
                kappa = kappaLoop;
                found_solution = true;
                break;
            end
        end

        kappaLoop = kappaLoop / kappaDecay;
    end

    if ~found_solution
        error('No feasible PROM solution was found.');
    end
end

function vec = replaceInfWithFinite(vec)
    if isempty(vec)
        return;
    end
    finiteVals = vec(isfinite(vec));
    if isempty(finiteVals)
        vec(:) = 0;
        return;
    end
    replacement = min(finiteVals);
    vec(isinf(vec)) = replacement;
    vec(isnan(vec)) = 0;
end

function out = toStringArray(x)
    if isempty(x)
        out = strings(0,1);
    else
        out = string(x(:));
    end
    out = out(strlength(out) > 0);
end
