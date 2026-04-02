function threshold_values = chooseThresholdValues(TF_exp, z_vector, thresholdSource)
    switch lower(string(thresholdSource))
        case "z_vector"
            threshold_values = unique(sort(z_vector(~isnan(z_vector) & isfinite(z_vector))));
            if isempty(threshold_values)
                threshold_values = 0;
            end
        otherwise
            threshold_values = unique(sort(TF_exp(~isnan(TF_exp) & isfinite(TF_exp))));
            if isempty(threshold_values)
                threshold_values = 0;
            end
    end
end