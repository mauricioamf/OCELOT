
function results = runKnockout(model, expression_data, beta_df, params)
% runOCELOTTFKnockout
% Self-contained OCELOT TF-knockout simulation workflow.
%
% This function runs TF knockout simulations and can optionally perform
% manuscript-style essentiality benchmarking after the KO simulations.
%
% Required inputs
%   model            COBRA model structure
%   expression_data  table with genes in RowNames and samples in VariableNames
%   beta_df          table with target genes in RowNames and TFs in VariableNames
%   params           struct with model, medium, and OCELOT settings
%
% Core required params
%   .selectedReference   WT reference condition present in expression_data
%   .biomassRxn          biomass reaction ID
%
% Growth-mode params
%   .growthType          'heterotrophic' (default), 'mixotrophic', or 'autotrophic'
%
% Heterotrophic / mixotrophic
%   .mediumRxns          exchange reactions opened in the medium
%   .glucoseRxn          glucose exchange reaction ID
%   .glucoseUptake       default -10
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
% Optional OCELOT params
%   .performEssentialityEvaluation default true
%   .essentialityTable   table with gene IDs and essentiality labels
%   .essentialityGeneColumn  optional gene-column name in essentialityTable
%   .essentialityLabelColumn optional label-column name in essentialityTable
%   .cutoffPerc         default 0.1
%   .showPlots          default 0; accepts 0/1/2/'cutoff'
%   .conditionMode       'selected' (default) or 'all'
%   .conditionList       explicit conditions to simulate
%   .thresholdSource     'tf_expression' or 'z_vector'
%   .thresholdSelection  'Growth', 'RMSE', 'Discriminability', or 'Composite'
%   .negativeBetaMode    'zero' (default) or 'keep'
%   .mu                  default 1
%   .maxInfeas           default 3
%   .kappaStart          default 1000
%   .kappaMin            default depends on growthType
%   .kappaDecay          default 10
%   .removeLastExchange  default true
%   .protectedRxns       reactions excluded from regulatory scaling
%   .additionalBoundChanges struct array with fields .rxn .value .bound
%
% Output
%   results  struct with WT threshold selection and TF knockout simulations

    params = applyDefaults(params);

    [avg_expr_data, unique_conditions, selectedWTReference] = preprocessExpressionData(expression_data, params.selectedReference);
    [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureBaseModel(model, params);
    reg = buildRegulatoryContext(modelWT_FBA, avg_expr_data, beta_df, selectedWTReference, params);

    wt = runWTThresholdSelection(modelWT_FBA, biomass_idx, reg, protectedRxnIdx, params);
    [best_threshold, best_idx, best_thr_value] = selectBestThreshold( ...
        wt.rmse_table, wt.threshold_names, wt.threshold_values, params.thresholdSelection);

    WT_biomass = wt.WT_results.(best_threshold).biomass;
    WT_kappa   = wt.WT_results.(best_threshold).kappa;

    if isfield(params, 'conditionList') && ~isempty(params.conditionList)
        condition_list = toCellstr(params.conditionList);
    elseif strcmpi(params.conditionMode, 'all')
        condition_list = unique_conditions;
    else
        condition_list = {selectedWTReference};
    end

    all_results = struct();
    evaluation = struct('performed', false, 'metrics_all', table(), 'metrics_cutoff', table(), ...
        'ocelot_eval_genes', {{}}, 'ocelot_y_true', [], 'covered_TF_list', {{}}, 'essentialityTableMatched', table());

    for c = 1:numel(condition_list)
        selectedCondition = condition_list{c};
        fprintf('\n=== Running condition: %s ===\n', selectedCondition);

        if ~ismember(selectedCondition, avg_expr_data.Properties.VariableNames)
            warning('Condition "%s" not found in averaged expression data. Skipping.', selectedCondition);
            continue;
        end

        condReg = buildConditionRegState(reg, avg_expr_data, beta_df, selectedCondition, best_thr_value, params);
        y_base = reg.gene_to_rxn_map * condReg.y_tgt_base;
        y_base = replaceInfWithFinite(y_base);

        yTable = array2table(y_base, 'RowNames', modelWT_FBA.rxns, 'VariableNames', {'WT'});
        for i = 1:numel(reg.TF_list)
            zko = condReg.z_bin;
            zko(i) = 0;
            yk = reg.gene_to_rxn_map * (reg.beta_matrix * (condReg.TF_exp .* zko));
            yk = replaceInfWithFinite(yk);
            yTable.(reg.TF_list{i}) = yk;
        end

        tauWT = max(abs(yTable{:,'WT'}));
        if tauWT == 0
            warning('WT tau is zero for condition %s. Skipping.', selectedCondition);
            continue;
        end

        y_KO_table = array2table(yTable{:,:} ./ tauWT, ...
            'RowNames', modelWT_FBA.rxns, ...
            'VariableNames', yTable.Properties.VariableNames);

        has_coverage = false(numel(reg.TF_list), 1);
        for i = 1:numel(reg.TF_list)
            delta = y_KO_table{:,'WT'} - y_KO_table{:, reg.TF_list{i}};
            has_coverage(i) = any(abs(delta) > 1e-9);
        end

        TF_growth_KO = table('Size', [numel(reg.TF_list), 1], ...
            'VariableTypes', {'double'}, ...
            'RowNames', reg.TF_list, ...
            'VariableNames', {'Biomass'});
        KO_fluxes = struct();

        for i = 1:numel(reg.TF_list)
            TF = reg.TF_list{i};
            [biomass, flux] = solveKOByRatio(modelWT_FBA, biomass_idx, WT_kappa, params.mu, ...
                reg.gene_to_rxn_map, y_KO_table, TF, protectedRxnIdx);

            TF_growth_KO{TF,'Biomass'} = biomass;
            if ~isempty(flux)
                KO_fluxes.(matlab.lang.makeValidName(TF)) = flux;
            end
        end

        condField = matlab.lang.makeValidName(selectedCondition);
        all_results.(condField).condition         = selectedCondition;
        all_results.(condField).KO_growth         = TF_growth_KO;
        all_results.(condField).KO_fluxes         = KO_fluxes;
        all_results.(condField).y_KO_table        = y_KO_table;
        all_results.(condField).has_coverage      = has_coverage;
        all_results.(condField).covered_TF_list   = reg.TF_list(has_coverage);
        all_results.(condField).WT_biomass        = WT_biomass;
        all_results.(condField).WT_kappa          = WT_kappa;
        all_results.(condField).threshold_value   = best_thr_value;
        all_results.(condField).threshold_name    = best_threshold;
        all_results.(condField).z_bin             = condReg.z_bin;
        all_results.(condField).TF_exp            = condReg.TF_exp;
    end

    if params.performEssentialityEvaluation
        evaluation = runEssentialityEvaluation(all_results, selectedWTReference, params);
    end

    results = struct();
    results.paramsUsed           = params;
    results.uniqueConditions     = unique_conditions;
    results.selectedWTReference  = selectedWTReference;
    results.avgExprData          = avg_expr_data;
    results.modelWT_FBA          = modelWT_FBA;
    results.FBAsol               = FBAsol;
    results.biomassWT_FBA        = biomassWT_FBA;
    results.protectedRxns        = protectedRxns;
    results.protectedRxnIdx      = protectedRxnIdx;
    results.WT_results           = wt.WT_results;
    results.rmse_table           = wt.rmse_table;
    results.best_threshold       = best_threshold;
    results.best_threshold_value = best_thr_value;
    results.best_threshold_index = best_idx;
    results.WT_kappa             = WT_kappa;
    results.regulatoryContext    = reg;
    results.all_results          = all_results;
    results.evaluation           = evaluation;
    results.metrics_all          = evaluation.metrics_all;
    results.metrics_cutoff       = evaluation.metrics_cutoff;
    results.covered_TF_list      = evaluation.covered_TF_list;
    results.ocelot_eval_genes    = evaluation.ocelot_eval_genes;
    results.ocelot_y_true        = evaluation.ocelot_y_true;
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
    defaults.conditionMode      = 'selected';
    defaults.negativeBetaMode   = 'zero';
    defaults.protectedRxns      = strings(0,1);
    defaults.additionalBoundChanges = struct([]);
    defaults.performEssentialityEvaluation = true;
    defaults.essentialityTable = [];
    defaults.essentialityGeneColumn = '';
    defaults.essentialityLabelColumn = '';
    defaults.cutoffPerc = 0.1;
    defaults.showPlots = 0;
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

    required = {'selectedReference','biomassRxn'};
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