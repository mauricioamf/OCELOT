function results = runEngineering(model, expression_data, beta_df, params)
% runEngineering
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
        params.thresholdSource = 'tf_expression';
    end

    if strlength(string(params.thresholdSelection)) == 0
        params.thresholdSelection = 'Composite';
    end

    if isempty(params.kappaMin)
        params.kappaMin = 1e-12;
    end

    if strlength(string(params.upsilonMode)) == 0
        params.upsilonMode = 'global_max';
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