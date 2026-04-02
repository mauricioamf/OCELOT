function condReg = buildConditionRegState(reg, avg_expr_data, beta_df, selectedCondition, thresholdValue, params)
    TF_exp = avg_expr_data{reg.TF_list, selectedCondition};
    target_exp_KO = avg_expr_data{beta_df.Properties.RowNames, selectedCondition};
    target_exp_valid_KO = target_exp_KO(reg.validTargets);

    beta_prime_KO = reg.beta_matrix .* TF_exp';
    z_vector_KO = pinv(beta_prime_KO) * target_exp_valid_KO;

    switch lower(string(params.thresholdSource))
        case "z_vector"
            z_bin = double(z_vector_KO >= thresholdValue);
        otherwise
            z_bin = double(TF_exp >= thresholdValue);
    end

    condReg = struct();
    condReg.TF_exp = TF_exp;
    condReg.z_vector = z_vector_KO;
    condReg.z_bin = z_bin;
    condReg.y_tgt_base = reg.beta_matrix * (TF_exp .* z_bin);
end