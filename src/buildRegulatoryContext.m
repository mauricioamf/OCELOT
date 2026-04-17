function reg = buildRegulatoryContext(modelWT_FBA, avg_expr_data, beta_df, selectedWTReference, params)
    TF_list = cellstr(string(beta_df.Properties.VariableNames(:)));
    target_genes = beta_df.Properties.RowNames;

    TF_exp = avg_expr_data{TF_list, selectedWTReference};
    target_exp = avg_expr_data{target_genes, selectedWTReference};  

    missingTFs = setdiff(TF_list, avg_expr_data.Properties.RowNames);
    if ~isempty(missingTFs)
        error('The following TFs are missing from expression_data.RowNames: %s', strjoin(missingTFs, ', '));
    end

    missingTargets = setdiff(target_genes, avg_expr_data.Properties.RowNames);
    if ~isempty(missingTargets)
        error('The following target genes are missing from expression_data.RowNames: %s', strjoin(missingTargets, ', '));
    end

    [~, targetGeneIdx_inModel] = ismember(target_genes, modelWT_FBA.genes);
    validTargets = targetGeneIdx_inModel > 0;

    targetGeneIdx_inModel = targetGeneIdx_inModel(validTargets);
    beta_df_valid = beta_df(validTargets, :);
    % target_gene_names_valid = target_genes(validTargets);
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