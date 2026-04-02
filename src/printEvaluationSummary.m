function printEvaluationSummary(metrics_cutoff)
    if isempty(metrics_cutoff)
        fprintf('\nNo evaluation results available at the selected cutoff.\n');
        return;
    end
    fprintf('\n=== Results restricted to metabolically covered TFs ===\n');
    fprintf('Mean coverage: %.1f%% of TFs affect at least one reaction\n', mean(metrics_cutoff.coverage_pct));
    fprintf('Mean essential (covered): %.1f | Mean non-essential (covered): %.1f\n', ...
        mean(metrics_cutoff.n_ess_covered), mean(metrics_cutoff.n_noness_covered));
    fprintf('\nAccuracy = %.3f | Sensitivity = %.3f | Specificity = %.3f | AUC = %.3f\n', ...
        mean(metrics_cutoff.Accuracy), mean(metrics_cutoff.Sensitivity), ...
        mean(metrics_cutoff.Specificity), mean(metrics_cutoff.AUC));
end