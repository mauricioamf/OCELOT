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