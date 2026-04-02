function evaluation = runEssentialityEvaluation(all_results, selectedWTReference, params)
    evaluation = struct('performed', false, 'metrics_all', table(), 'metrics_cutoff', table(), ...
        'ocelot_eval_genes', {{}}, 'ocelot_y_true', [], 'covered_TF_list', {{}}, 'essentialityTableMatched', table());

    if ~isfield(params, 'essentialityTable') || isempty(params.essentialityTable)
        warning('performEssentialityEvaluation is true but params.essentialityTable is missing or empty. Skipping evaluation.');
        return;
    end

    essentialityTable = params.essentialityTable;
    [geneCol, labelCol] = resolveEssentialityColumns(essentialityTable, params);
    cutoff_range = 0.05:0.05:1.0;
    metrics_all = table();
    selected_eval_genes = {};
    selected_y_true = [];
    selected_covered_tfs = {};
    selected_matched = table();

    condFields = fieldnames(all_results);
    for i = 1:numel(condFields)
        condField = condFields{i};
        R = all_results.(condField);
        if ~isfield(R, 'has_coverage') || ~isfield(R, 'KO_growth')
            continue;
        end

        covered_TF_list = R.covered_TF_list;
        if isempty(covered_TF_list)
            fprintf('  No covered TFs for condition %s. Skipping essentiality evaluation.\n', R.condition);
            continue;
        end

        [lia, loc] = ismember(string(essentialityTable.(geneCol)), string(covered_TF_list));
        if sum(lia) == 0
            fprintf('  No covered TFs match essentiality table for condition %s.\n', R.condition);
            continue;
        end

        y_true = essentialityTable.(labelCol)(lia);
        if islogical(y_true)
            y_true = double(y_true);
        end
        biomass_vals = R.KO_growth{covered_TF_list, 'Biomass'};
        biomass_vals = biomass_vals(loc(lia));
        eval_genes = cellstr(string(essentialityTable.(geneCol)(lia)));
        matchedTable = essentialityTable(lia, :);

        if all(y_true == 0) || all(y_true == 1)
            fprintf('  Only one class present in essentiality labels for condition %s. Skipping.\n', R.condition);
            continue;
        end

        n_ess_covered = sum(y_true == 1);
        n_noness_covered = sum(y_true == 0);
        fprintf('  Covered TFs in essentiality table for %s: %d essential, %d non-essential\n', ...
            R.condition, n_ess_covered, n_noness_covered);

        scores_continuous = 1 - (biomass_vals ./ R.WT_biomass);
        [fpRate, tpRate, ~, AUC] = perfcurve(y_true, scores_continuous, 1);

        for cutoffVal = cutoff_range
            cutoff = cutoffVal * R.WT_biomass;
            pred_viable = biomass_vals >= cutoff;
            y_pred_ess = double(~pred_viable);

            C = confusionmat(y_true, y_pred_ess);
            if ~isequal(size(C), [2 2])
                continue;
            end
            TN = C(1,1); FP = C(1,2);
            FN = C(2,1); TP = C(2,2);

            accuracy = (TP + TN) / sum(C(:));
            sensitivity = TP / (TP + FN + eps);
            specificity = TN / (TN + FP + eps);

            metrics_loop = table(string(R.condition), cutoffVal, accuracy, sensitivity, specificity, AUC, ...
                'VariableNames', {'Condition','cutoffVal','Accuracy','Sensitivity','Specificity','AUC'});
            metrics_loop.tpRate{1} = tpRate;
            metrics_loop.fpRate{1} = fpRate;
            metrics_loop.C{1} = C;
            metrics_loop.n_ess_covered(1) = n_ess_covered;
            metrics_loop.n_noness_covered(1) = n_noness_covered;
            metrics_loop.coverage_pct(1) = sum(R.has_coverage) / numel(R.has_coverage) * 100;
            metrics_all = [metrics_all; metrics_loop]; %#ok<AGROW>
        end

        if strcmp(R.condition, selectedWTReference)
            selected_eval_genes = eval_genes;
            selected_y_true = y_true;
            selected_covered_tfs = covered_TF_list;
            selected_matched = matchedTable;
        end
    end

    metrics_all = rmmissing(metrics_all);
    metrics_cutoff = table();
    if ~isempty(metrics_all)
        metrics_cutoff = metrics_all(metrics_all.cutoffVal == params.cutoffPerc, :);
    end

    if isempty(selected_eval_genes) && ~isempty(condFields)
        for i = 1:numel(condFields)
            condField = condFields{i};
            R = all_results.(condField);
            covered_TF_list = R.covered_TF_list;
            [lia, ~] = ismember(string(essentialityTable.(geneCol)), string(covered_TF_list));
            if sum(lia) > 0
                selected_eval_genes = cellstr(string(essentialityTable.(geneCol)(lia)));
                selected_y_true = essentialityTable.(labelCol)(lia);
                selected_covered_tfs = covered_TF_list;
                selected_matched = essentialityTable(lia, :);
                break;
            end
        end
    end

    evaluation.performed = ~isempty(metrics_all);
    evaluation.metrics_all = metrics_all;
    evaluation.metrics_cutoff = metrics_cutoff;
    evaluation.ocelot_eval_genes = selected_eval_genes;
    evaluation.ocelot_y_true = selected_y_true;
    evaluation.covered_TF_list = selected_covered_tfs;
    evaluation.essentialityTableMatched = selected_matched;

    if evaluation.performed
        printEvaluationSummary(metrics_cutoff);
        plotEvaluation(metrics_cutoff, selectedWTReference, params.showPlots);
    end
end