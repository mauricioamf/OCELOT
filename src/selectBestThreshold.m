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