function results = OCELOT(model, expression_data, beta_df, params)
% OCELOT
% Unified OCELOT workflow for TF knockout simulations and engineering.
%
% This wrapper unifies the TF-knockout and engineering workflows into a
% single entry point, controlled by flags in params.
%
% Required inputs
%   model            COBRA model structure
%   expression_data  table with genes in RowNames and samples in VariableNames
%   beta_df          table with target genes in RowNames and TFs in VariableNames
%   params           struct with OCELOT settings
%
% Task flags (both optional)
%   .runKnockout   default true
%   .runEngineering  default false
%
% Evaluation / validation flags
%   .performEssentialityEvaluation   default true  (used by KO workflow)
%   .validateStrategies              default true  (used by engineering workflow)
%
% Notes
% - This unified function intentionally reuses the dedicated OCELOT task
%   functions so the existing KO evaluations and engineering validations are
%   preserved without changing their internal implementations.
% - For TF knockout only, engineering-only parameters such as .prodRxn are
%   not required.
% - For engineering, supply the same parameters expected by
%   runOCELOTEngineering.
%
% Output
%   results struct with fields:
%       .paramsUsed
%       .tasks
%       .knockout     result from runOCELOTTFKnockout (if requested)
%       .engineering  result from runOCELOTEngineering (if requested)
%
% Example
%   params.runKnockout = true;
%   params.runEngineering = true;
%   params.performEssentialityEvaluation = true;
%   params.validateStrategies = true;
%   results = runOCELOTUnified(model, expression_data, beta_df, params);

    if nargin < 4 || isempty(params)
        params = struct();
    end

    params = applyUnifiedDefaults(params);

    if ~params.runKnockout && ~params.runEngineering
        error('At least one task must be enabled: params.runKnockout or params.runEngineering.');
    end

    results = struct();
    results.paramsUsed = params;
    results.tasks = struct( ...
        'runKnockout', params.runKnockout, ...
        'runEngineering', params.runEngineering, ...
        'performEssentialityEvaluation', params.performEssentialityEvaluation, ...
        'validateStrategies', params.validateStrategies);

    % Run TF knockout workflow (including optional manuscript-style evaluation)
    if params.runKnockout
        koParams = params;
        koParams.performEssentialityEvaluation = params.performEssentialityEvaluation;
        results.knockout = runKnockout(model, expression_data, beta_df, koParams);
    else
        results.knockout = struct();
    end

    % Run engineering workflow (including optional validation)
    if params.runEngineering
        engParams = params;
        engParams.validateStrategies = params.validateStrategies;
        results.engineering = runEngineering(model, expression_data, beta_df, engParams);
    else
        results.engineering = struct();
    end

    % Convenience aliases to reduce downstream code changes
    if params.runKnockout
        results.WT_results_KO         = results.knockout.WT_results;
        results.rmse_table_KO         = results.knockout.rmse_table;
        results.best_threshold_KO     = results.knockout.best_threshold;
        results.best_threshold_value_KO = results.knockout.best_threshold_value;
        results.all_results_KO        = results.knockout.all_results;
        results.metrics_all           = results.knockout.metrics_all;
        results.metrics_cutoff        = results.knockout.metrics_cutoff;
        results.covered_TF_list       = results.knockout.covered_TF_list;
        results.ocelot_eval_genes     = results.knockout.ocelot_eval_genes;
        results.ocelot_y_true         = results.knockout.ocelot_y_true;
    end

    if params.runEngineering
        results.WT_results_ENG           = results.engineering.WT_results;
        results.rmse_table_ENG           = results.engineering.rmse_table;
        results.best_threshold_ENG       = results.engineering.best_threshold;
        results.best_threshold_value_ENG = results.engineering.best_threshold_value;
        results.resultsTable             = results.engineering.resultsTable;
        results.filteredResultsTable     = results.engineering.filteredResultsTable;
        results.validatedResultsTable    = results.engineering.validatedResultsTable;
    end
end

function params = applyUnifiedDefaults(params)
    defaults.runKnockout = true;
    defaults.runEngineering = false;
    defaults.performEssentialityEvaluation = true;
    defaults.validateStrategies = true;

    fn = fieldnames(defaults);
    for i = 1:numel(fn)
        if ~isfield(params, fn{i}) || isempty(params.(fn{i}))
            params.(fn{i}) = defaults.(fn{i});
        end
    end
end
