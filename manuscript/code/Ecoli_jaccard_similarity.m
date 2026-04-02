%% Specificity analysis across metabolite targets

% Choose target product reactions to compare
prod_rxn_list = {'EX_ac_e','EX_etoh_e','EX_for_e','EX_lac__D_e','EX_succ_e'};

nProd = numel(prod_rxn_list);

prod_idx_list  = zeros(nProd,1);
prod_name_list = strings(nProd,1);

for p = 1:nProd
    prod_idx_list(p) = find(strcmp(modelWT_FBA.rxns, prod_rxn_list{p}), 1);
    if isempty(prod_idx_list(p))
        error('Product reaction %s not found.', prod_rxn_list{p});
    end

    met_idx = find(model.S(:,prod_idx_list(p)));
    if ~isempty(met_idx)
        prod_name_list(p) = string(model.metNames{met_idx(1)});
    else
        prod_name_list(p) = string(prod_rxn_list{p});
    end
end

% Build TF -> affected reaction sets
tf_rxn_sets = cell(nTF,1);
tf_reg_strength = zeros(nRxns, nTF);

for j = 1:nTF
    affected_idx = find(abs(upsilon(:,j)) > 1e-12);
    tf_rxn_sets{j} = affected_idx;
    tf_reg_strength(:,j) = abs(upsilon(:,j));
end

% Pairwise Jaccard overlap between TF-affected reaction sets
tf_jaccard = NaN(nTF, nTF);

for i = 1:nTF
    A = tf_rxn_sets{i};
    for j = 1:nTF
        B = tf_rxn_sets{j};

        if isempty(A) && isempty(B)
            tf_jaccard(i,j) = 1;
        else
            inter = numel(intersect(A,B));
            uni   = numel(union(A,B));
            tf_jaccard(i,j) = inter / max(uni,1);
        end
    end
end

tf_jaccard_table = array2table(tf_jaccard, ...
    'VariableNames', TF_list, ...
    'RowNames', TF_list);

% Optional distance form
tf_jaccard_distance = 1 - tf_jaccard;

% Solve FBA for each target metabolite
active_rxn_mat           = false(nRxns, nProd);
regulated_active_rxn_mat = false(nRxns, nProd);
flux_mat                 = zeros(nRxns, nProd);

for p = 1:nProd
    p_idx = prod_idx_list(p);

    modelP = modelWT_FBA;
    modelP.c(:) = 0;
    modelP.c(p_idx) = 1;
    modelP.lb(p_idx) = -1000;
    modelP.ub(p_idx) = 1000;
    modelP.lb(biomass_idx) = biomassWT_FBA * 0.5;

    solP = optimizeCbModel(modelP, 'max');

    if isempty(solP.f) || isnan(solP.f)
        warning('No feasible solution for %s', prod_rxn_list{p});
        continue;
    end

    flux_mat(:,p) = solP.x;

    active_idx = find(abs(solP.x) > 1e-8);
    reg_active_idx = active_idx(any(abs(upsilon(active_idx,:)) > 1e-12, 2));

    active_rxn_mat(active_idx,p) = true;
    regulated_active_rxn_mat(reg_active_idx,p) = true;
end

% Pairwise Jaccard overlaps between metabolites
met_jaccard_active = NaN(nProd, nProd);
met_jaccard_reg    = NaN(nProd, nProd);

for i = 1:nProd
    Ai  = find(active_rxn_mat(:,i));
    Ari = find(regulated_active_rxn_mat(:,i));

    for j = 1:nProd
        Aj  = find(active_rxn_mat(:,j));
        Arj = find(regulated_active_rxn_mat(:,j));

        if isempty(Ai) && isempty(Aj)
            met_jaccard_active(i,j) = 1;
        else
            met_jaccard_active(i,j) = numel(intersect(Ai,Aj)) / max(numel(union(Ai,Aj)),1);
        end

        if isempty(Ari) && isempty(Arj)
            met_jaccard_reg(i,j) = 1;
        else
            met_jaccard_reg(i,j) = numel(intersect(Ari,Arj)) / max(numel(union(Ari,Arj)),1);
        end
    end
end

met_active_table = array2table(met_jaccard_active, ...
    'VariableNames', cellstr(prod_name_list), ...
    'RowNames', cellstr(prod_name_list));

met_reg_table = array2table(met_jaccard_reg, ...
    'VariableNames', cellstr(prod_name_list), ...
    'RowNames', cellstr(prod_name_list));

% TF coverage of each metabolite route
tf_route_coverage = zeros(nTF, nProd);
tf_route_weighted = zeros(nTF, nProd);

for p = 1:nProd
    route_idx = find(active_rxn_mat(:,p));
    if isempty(route_idx)
        continue;
    end

    % Unweighted coverage: fraction of active reactions touched by TF
    for j = 1:nTF
        tf_route_coverage(j,p) = sum(abs(upsilon(route_idx,j)) > 1e-12) / numel(route_idx);
    end

    % Weighted coverage: absolute flux-weighted regulatory coverage
    flux_weights = abs(flux_mat(route_idx,p));
    flux_weights = flux_weights / (sum(flux_weights) + eps);

    for j = 1:nTF
        tf_route_weighted(j,p) = sum((abs(upsilon(route_idx,j)) > 1e-12) .* flux_weights);
    end
end

tf_route_coverage_table = array2table(tf_route_coverage, ...
    'VariableNames', cellstr(prod_name_list), ...
    'RowNames', TF_list);

tf_route_weighted_table = array2table(tf_route_weighted, ...
    'VariableNames', cellstr(prod_name_list), ...
    'RowNames', TF_list);

% Shared vs unique active/regulated reactions
shared_active_all = all(active_rxn_mat, 2);
shared_reg_all    = all(regulated_active_rxn_mat, 2);

unique_active_counts = zeros(nProd,1);
unique_reg_counts    = zeros(nProd,1);
shared_active_counts = zeros(nProd,1);
shared_reg_counts    = zeros(nProd,1);

for p = 1:nProd
    this_active = active_rxn_mat(:,p);
    this_reg    = regulated_active_rxn_mat(:,p);

    other_active = any(active_rxn_mat(:, setdiff(1:nProd,p)), 2);
    other_reg    = any(regulated_active_rxn_mat(:, setdiff(1:nProd,p)), 2);

    unique_active_counts(p) = sum(this_active & ~other_active);
    unique_reg_counts(p)    = sum(this_reg & ~other_reg);
    shared_active_counts(p) = sum(this_active & shared_active_all);
    shared_reg_counts(p)    = sum(this_reg & shared_reg_all);
end

route_overlap_summary = table( ...
    prod_name_list, ...
    sum(active_rxn_mat,1)', ...
    sum(regulated_active_rxn_mat,1)', ...
    shared_active_counts, ...
    shared_reg_counts, ...
    unique_active_counts, ...
    unique_reg_counts, ...
    'VariableNames', {'Metabolite','ActiveRxns','RegulatedActiveRxns', ...
                      'SharedActiveRxns','SharedRegulatedRxns', ...
                      'UniqueActiveRxns','UniqueRegulatedRxns'});

% Reaction-level summary for each metabolite
reaction_summary = table();

for p = 1:nProd
    active_idx = find(active_rxn_mat(:,p));
    if isempty(active_idx)
        continue;
    end

    tmp = table( ...
        repmat(prod_name_list(p), numel(active_idx), 1), ...
        modelWT_FBA.rxns(active_idx), ...
        string(modelWT_FBA.rxnNames(active_idx)), ...
        flux_mat(active_idx,p), ...
        regulated_active_rxn_mat(active_idx,p), ...
        shared_active_all(active_idx), ...
        shared_reg_all(active_idx), ...
        'VariableNames', {'Metabolite','ReactionID','ReactionName','Flux', ...
                          'IsRegulated','SharedAcrossAllProducts','SharedRegulatedAcrossAllProducts'});

    reaction_summary = [reaction_summary; tmp];
end

% Top TFs per metabolite route
top_n = min(10, nTF);
top_tf_per_met = table();

for p = 1:nProd
    [sorted_cov, ord] = sort(tf_route_weighted(:,p), 'descend');
    keep = ord(1:top_n);

    tmp = table( ...
        repmat(prod_name_list(p), top_n, 1), ...
        string(TF_list(keep)), ...
        tf_route_coverage(keep,p), ...
        tf_route_weighted(keep,p), ...
        'VariableNames', {'Metabolite','TF','CoverageFraction','FluxWeightedCoverage'});

    top_tf_per_met = [top_tf_per_met; tmp];
end

% Jaccard overlap among TFs within predicted KO set(s)
if height(resultsTableFilter) > 0
    % Collect all unique TFs used in predicted strategies
    all_predicted_tfs = {};
    for i = 1:height(resultsTableFilter)
        if ~isempty(resultsTableFilter.Knockouts{i})
            all_predicted_tfs = [all_predicted_tfs; resultsTableFilter.Knockouts{i}(:)];
        end
        if ~isempty(resultsTableFilter.Activated{i})
            all_predicted_tfs = [all_predicted_tfs; resultsTableFilter.Activated{i}(:)];
        end
    end
    all_predicted_tfs = unique(all_predicted_tfs);

    pred_idx = find(ismember(TF_list, all_predicted_tfs));
    pred_tf_names = TF_list(pred_idx);

    pred_tf_jaccard = tf_jaccard(pred_idx, pred_idx);
    pred_tf_jaccard_table = array2table(pred_tf_jaccard, ...
        'VariableNames', pred_tf_names, ...
        'RowNames', pred_tf_names);
else
    pred_tf_jaccard_table = table();
end

% Export
writetable(route_overlap_summary, 'Ecoli_route_overlap_summary.xlsx');
writetable(tf_route_coverage_table, 'Ecoli_tf_route_coverage.xlsx', 'WriteRowNames', true);
writetable(tf_route_weighted_table, 'Ecoli_tf_route_weighted_coverage.xlsx', 'WriteRowNames', true);
writetable(met_active_table, 'Ecoli_metabolite_active_jaccard.xlsx', 'WriteRowNames', true);
writetable(met_reg_table, 'Ecoli_metabolite_regulated_jaccard.xlsx', 'WriteRowNames', true);
writetable(tf_jaccard_table, 'Ecoli_tf_reaction_jaccard.xlsx', 'WriteRowNames', true);
writetable(reaction_summary, 'Ecoli_reaction_summary_by_metabolite.xlsx');
writetable(top_tf_per_met, 'Ecoli_topTFs_per_metabolite.xlsx');
if ~isempty(pred_tf_jaccard_table)
    writetable(pred_tf_jaccard_table, 'Ecoli_predictedTF_jaccard.xlsx', 'WriteRowNames', true);
end