%% Configurations for running on the HPC
% Initialize parallel pool for HPC
poolobj = gcp('nocreate');
if isempty(poolobj)
   % Use SLURM CPU count if available, otherwise default to maximum local threads
   numWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(numWorkers), numWorkers = maxNumCompThreads; end
   parpool('local', numWorkers);
end

initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
pctRunOnAll initCobraToolbox(false);
pctRunOnAll changeCobraSolver('gurobi','LP',0);
changeCobraSolverParams('LP', 'feasTol', 1e-9);
clc;clear
tic

%% Load metabolic model and other data

% Load metabolic model
load('../models/yeast-GEM.mat');

% Load expression data
expression_data = readtable('../data/Sce_count_matrix.csv', 'ReadRowNames', true, 'VariableNamingRule', 'preserve');

% Load beta matrix as dataframe
beta_df = readtable('../data/Sce_GRN_lasso.csv', 'ReadRowNames', true);

% Load gene essentiality
essentialityTable = readtable('../data/Sce_gene_essentiality.csv');

%% Set parameters
warning("on");
selectedReference = "FY4_6"; % Cluster_1 or control__wt_glc or MG1655
showPlots = "cutoff"; % 0 = no plots, 1 = confusion matrix and ROC only, 2 = all
glucoseUptake = -10;
mu = 1;
thresholdSelection = 'RMSE'; % 'RMSE' or 'Growth'
biomass_rxn = 'r_2111';
cutoffPerc = 0.1;
max_infeas   = 3;  % stop after n consecutive infeasible thresholds

%% Process expression data
all_samples = expression_data.Properties.VariableNames;
condition_names = regexprep(all_samples, '__\d+$', '');  % remove __1, __2, etc.

for i = 1:numel(condition_names)
    if ~isempty(condition_names{i}) && isstrprop(condition_names{i}(1), 'digit')
        condition_names{i} = ['X' condition_names{i}];
    end
end

unique_conditions = unique(condition_names, 'stable');
wt_idx = find(contains(unique_conditions, selectedReference), 1, 'first');
if isempty(wt_idx), error('WT reference not found in unique_conditions'); end
selectedWTReference = unique_conditions{wt_idx};

avg_expr_data = array2table( ...
    zeros(size(expression_data,1), numel(unique_conditions)), ...
    'RowNames', expression_data.Properties.RowNames, ...
    'VariableNames', unique_conditions);

for c = 1:numel(unique_conditions)
    cond = unique_conditions{c};
    idx  = strcmp(condition_names, cond);  % all replicates for this condition
    avg_expr_data{:,c} = mean(expression_data{:,idx}, 2, 'omitnan');
end

%% Optimizations for WT
% FBA for WT growth rate
modelWT_FBA = model;

exchangeIndex = find(contains(modelWT_FBA.rxnNames, "exchange"));
exchangeIDs = modelWT_FBA.rxns(exchangeIndex);
exchangeIDs(end) = [];
biomass_idx = find(strcmp(modelWT_FBA.rxns, biomass_rxn));

M9_components = ["r_1654","r_1992", "r_2005","r_2060", ...
                 "r_1861","r_1832","r_2100","r_4593", ...
                 "r_4595","r_4596","r_4597","r_2049", ...
                 "r_4594","r_4600","r_2020"];

modelWT_FBA = changeRxnBounds(modelWT_FBA, exchangeIDs, 0, 'l');
modelWT_FBA = changeRxnBounds(modelWT_FBA, M9_components, -1000, 'l');
modelWT_FBA = changeObjective(modelWT_FBA, modelWT_FBA.rxns{biomass_idx});

modelWT_FBA = changeRxnBounds(modelWT_FBA, 'r_1714', glucoseUptake, 'b'); %D-glucose
FBAsol = optimizeCbModel(modelWT_FBA);
biomassWT_FBA = FBAsol.f; % Store the FBA wild-type biomass value

%% Find WT kappa and expression threshold
% Calculate WT parameters and for multiple z thresholds
TF_list = beta_df.Properties.VariableNames';  % T x 1

TF_exp = avg_expr_data{TF_list, selectedWTReference};     % returns numeric column vector
target_exp = avg_expr_data{beta_df.Properties.RowNames, selectedWTReference};

assert(isequal(size(TF_exp), [numel(TF_list),1]), 'TF_exp must be T×1');
assert(isequal(size(target_exp), [height(beta_df),1]), 'target_exp must be G×1');

[n_targets, n_TFs] = size(beta_df);
target_genes = beta_df.Properties.RowNames;

% Map target genes to model genes
[~, targetGeneIdx_inModel] = ismember(target_genes, modelWT_FBA.genes);
validTargets = targetGeneIdx_inModel > 0;

targetGeneIdx_inModel = targetGeneIdx_inModel(validTargets);
beta_df_valid = beta_df(validTargets, :);
target_gene_names_valid = target_genes(validTargets);
target_exp_valid = target_exp(validTargets);

modelWT_FBA = buildRxnGeneMat(modelWT_FBA);
gene_to_rxn_map = modelWT_FBA.rxnGeneMat(:, targetGeneIdx_inModel);  % (n_rxns x n_validTargets)

% Build matrices and vectors
beta_matrix = beta_df_valid{:,:};
beta_matrix(beta_matrix < 0) = 0;
beta_prime = beta_matrix .* TF_exp';

z_vector = pinv(beta_prime) * target_exp_valid;

% Thresholds to test
thresholds = unique(sort(z_vector));   % one threshold per TF expression level

% Output structs
y_all_struct   = struct();  % holds y_all_tau per threshold
z_struct       = struct();  % holds z_vec per threshold
z_struct_cont  = struct();  % holds z_vec per threshold

% Loop through thresholds
for t = 1:numel(thresholds)
    thr     = thresholds(t);
    z_bin   = double(z_vector >= thr);
    z_cont  = z_vector;
    z_cont(z_cont <= thr) = 0;
    [nG,nT] = size(beta_df);

    name  = sprintf('thr_expr_%g',thr);
    name  = matlab.lang.makeValidName(name);

    % WT reaction scalings
    TF_activity = TF_exp .* z_bin;
    y_tgt = beta_matrix * TF_activity;
    y_rxn = gene_to_rxn_map  * y_tgt;
    y_rxn(isinf(y_rxn)) = min(y_rxn(~isinf(y_rxn)));

    T = array2table(y_rxn, ...
        'RowNames',modelWT_FBA.rxns,'VariableNames',{'WT'});

    % normalize
    tau = max(abs(T{:,'WT'}));
    T   = array2table(T{:,:}/tau, ...
        'RowNames',modelWT_FBA.rxns,'VariableNames',T.Properties.VariableNames);

    y_all_struct.(name)  = T;
    z_struct.(name)      = z_bin;
    z_struct_cont.(name) = z_cont;
end

TF_growth_WT = table('Size', [length(TF_list)+1 1], ...
    'VariableTypes', {'double'}, ...
    'RowNames', ['WT'; TF_list], ...
    'VariableNames', {'Biomass'});

% Get WT solution across thresholds
WT_results = struct();
threshold_names = fieldnames(y_all_struct);
infeas_count = 0;  % counter for consecutive infeasible runs

for tr = 1:length(threshold_names)
    threshold_name = threshold_names{tr};

    fprintf('\n\n[%s] --- WT PROM for threshold: %s ---', ...
        datetime('now','Format','HH:mm:ss'), ...
        threshold_name);

    try
        modelWT = modelWT_FBA;
        S_WT = sparse(modelWT.S);
        [m_WT, n_WT] = size(S_WT);
        lb_base_WT = modelWT.lb;
        ub_base_WT = modelWT.ub;

        TF_exp_vec = TF_exp;

        y_all_loop = y_all_struct.(threshold_name);
        yWT = y_all_loop{:,'WT'};
        z_bin = z_struct.(threshold_name);
        z_cont = z_struct_cont.(threshold_name);
        target_exp_hat = beta_prime * z_cont;

        valid = ~isnan(yWT);
        rxnIDs = y_all_loop.Properties.RowNames(valid);
        yWT_valid = yWT(valid);
        yWT_valid(~isfinite(yWT_valid)) = 1;

        found_solution = false;
        lbp_WT = lb_base_WT;
        ubp_WT = ub_base_WT;

        for j_WT = 1:numel(rxnIDs)
            rxn_WT = rxnIDs{j_WT};
            idx_WT = find(strcmp(modelWT.rxns, rxn_WT), 1);

            if isempty(idx_WT), continue; end
            if gene_to_rxn_map(idx_WT,:) == 0, continue; end % reaction is unregulated

            y_WT = yWT_valid(j_WT);

            if abs(lb_base_WT(idx_WT)) > 1e-6 && abs(ub_base_WT(idx_WT)) > 1e-6
                lbp_WT(idx_WT) = max(lb_base_WT(idx_WT) * y_WT, -1000);
                ubp_WT(idx_WT) = min(ub_base_WT(idx_WT) * y_WT, 1000);
            elseif lb_base_WT(idx_WT) < -1e-6
                lbp_WT(idx_WT) = max(lb_base_WT(idx_WT) * y_WT, -1000);
            elseif ub_base_WT(idx_WT) > 1e-6
                ubp_WT(idx_WT) = min(ub_base_WT(idx_WT) * y_WT, 1000);
            end
        end

        lbp_WT(~isfinite(lbp_WT)) = lb_base_WT(~isfinite(lbp_WT));
        ubp_WT(~isfinite(ubp_WT)) = ub_base_WT(~isfinite(ubp_WT));

        infeasible_idx = find(lbp_WT > ubp_WT);
        if ~isempty(infeasible_idx)
            fprintf('\n[%s] Found %d reactions with lb > ub after scaling — fixing them.', ...
                datetime('now','Format','HH:mm:ss'), numel(infeasible_idx));
            for k_inf = infeasible_idx'
                lbp_WT(k_inf) = min(lbp_WT(k_inf), ubp_WT(k_inf) - 1e-6);
                ubp_WT(k_inf) = max(lbp_WT(k_inf), ubp_WT(k_inf) + 1e-6);
            end
        end

        % Tune kappa
        kappa_WT = 1000;
        kappa_WT_min = 1e-12;

        while kappa_WT >= kappa_WT_min
            gp_WT.A = sparse([
                S_WT,             sparse(m_WT,n_WT),    sparse(m_WT,n_WT);
               -speye(n_WT),     -speye(n_WT),          sparse(n_WT,n_WT);
                speye(n_WT),      sparse(n_WT,n_WT),   -speye(n_WT)
            ]);
            gp_WT.rhs        = [ zeros(m_WT,1); -lbp_WT; ubp_WT ];
            gp_WT.sense      = [ repmat('=', m_WT, 1); repmat('<', n_WT, 1); repmat('<', n_WT, 1) ];
            gp_WT.modelsense = 'min';

            gp_WT.obj = [ zeros(n_WT,1); kappa_WT*ones(n_WT,1); kappa_WT*ones(n_WT,1) ];
            gp_WT.obj(biomass_idx) = -mu;

            gp_WT.lb = [lb_base_WT; zeros(2*n_WT,1)];
            gp_WT.ub = [ub_base_WT; inf(2*n_WT,1)];
            gp_WT.varnames = [strcat('v_', modelWT.rxns); strcat('α_', modelWT.rxns); strcat('β_', modelWT.rxns)];

            params = struct();
            params.OutputFlag = 0;
            result_WT = gurobi(gp_WT, params);

            if isfield(result_WT, 'x')
                alpha_vals = result_WT.x(n_WT+1 : 2*n_WT);
                beta_vals  = result_WT.x(2*n_WT+1 : 3*n_WT);
                num_slacks_used = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

                v_sol = result_WT.x(1:n_WT);
                biomass = v_sol(biomass_idx);

                if biomass > 1e-6
                    fprintf('\n[%s] → predicted biomass = %.4f, slacks = %d, kappa = %.3e', ...
                        datetime('now','Format','HH:mm:ss'), ...
                        biomass, num_slacks_used, kappa_WT);

                    R = struct();
                    R.biomass = biomass;
                    R.model = gp_WT;
                    R.result = result_WT;
                    R.kappa = kappa_WT;
                    found_solution = true;
                    break;
                end
            else
                fprintf('\n[%s] %s → No solution at kappa = %.3e', ...
                    datetime('now','Format','HH:mm:ss'), ...
                    threshold_name, kappa_WT);
            end
            kappa_WT = kappa_WT / 10;
        end

        if ~found_solution
            error('No feasible solution found for %s.', threshold_name);
        end

        R.target_expr = array2table(target_exp_hat, ...
            'VariableNames', {threshold_name}, ...
            'RowNames', beta_df_valid.Properties.RowNames);
        WT_results.(threshold_name) = R;

    catch
        infeas_count = infeas_count + 1;
        warning('No feasible solution found for %s.', threshold_name);

        R.biomass = NaN;
        R.model = [];
        R.result = [];
        R.kappa = NaN;

        R.target_expr = array2table(target_exp_hat, ...
            'VariableNames', {threshold_name}, ...
            'RowNames', beta_df_valid.Properties.RowNames);
        WT_results.(threshold_name) = R;

        if infeas_count >= max_infeas
            fprintf('\n[%s] Stopping WT simulations after %d consecutive infeasible thresholds.\n', ...
                datetime('now','Format','HH:mm:ss'), max_infeas);
            break;
        end
    end
end

%% Compute RMSE for calculated target_exp_hat
simulatedThresholds = fieldnames(WT_results);

rmse_table = table('Size', [length(simulatedThresholds), 4], ...
    'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Threshold', 'RMSE', 'R2', 'Growth'});

for th = 1:length(simulatedThresholds)
    threshold_name = simulatedThresholds{th};
    predicted_tbl = WT_results.(threshold_name).target_expr;
    predicted_vec = predicted_tbl{:,1};
    target_exp_valid = target_exp(validTargets);

    predicted_vec = log(abs(predicted_vec));
    target_exp_valid = log(abs(target_exp_valid));

    residual = target_exp_valid - predicted_vec;

    goodIdx = isfinite(residual);
    target_exp_valid = target_exp_valid(goodIdx);
    predicted_vec = predicted_vec(goodIdx);

    residual = target_exp_valid - predicted_vec;
    rmse_value = sqrt(mean(residual.^2));

    ss_res = sum((predicted_vec - target_exp_valid).^2);
    ss_tot = sum((target_exp_valid - mean(target_exp_valid)).^2);
    r2 = 1 - (ss_res / ss_tot);

    rmse_table.Threshold(th) = threshold_name;
    rmse_table.RMSE(th) = rmse_value;
    rmse_table.R2(th) = r2;
    rmse_table.Growth(th) = WT_results.(threshold_name).biomass;
end

valid_idx = ~isnan(rmse_table.RMSE) & ~isnan(rmse_table.Growth) & rmse_table.Growth > 0;

switch thresholdSelection
    case 'Growth'
        [~, best_local_idx] = max(rmse_table.Growth(valid_idx));
    case 'RMSE'
        [~, best_local_idx] = min(rmse_table.RMSE(valid_idx));
end

best_idx_all = find(valid_idx);
best_idx = best_idx_all(best_local_idx);
best_threshold = rmse_table.Threshold(best_idx);

if ~isempty(best_threshold)
    fprintf('\n[%s] Selected threshold: %s     Criteria: %s \n', ...
        datetime('now','Format','HH:mm:ss'), best_threshold, thresholdSelection);
else
    error('\nNo threshold could be selected: either RMSE and/or growth has NaN or no growth is higher than 0.');
end

%% Loop over all conditions

WT_biomass    = WT_results.(best_threshold).biomass;
WT_kappa      = WT_results.(best_threshold).kappa;

all_results = struct();
metrics_all = table();

for c = 1:numel(unique_conditions)
    selectedCondition = unique_conditions{c};

    fprintf('\n[%s] === Running condition: %s ===\n', ...
        datetime('now','Format','HH:mm:ss'), selectedCondition);

    % Calculate KO parameters
    TF_exp_KO = avg_expr_data{TF_list, selectedCondition};
    target_exp_KO = avg_expr_data{beta_df.Properties.RowNames, selectedCondition};

    assert(isequal(size(TF_exp_KO), [numel(TF_list),1]), 'TF_exp must be T×1');
    assert(isequal(size(target_exp_KO), [height(beta_df),1]), 'target_exp must be G×1');

    target_exp_valid_KO = target_exp_KO(validTargets);

    beta_matrix_KO = beta_matrix;
    beta_prime_KO = beta_matrix_KO .* TF_exp_KO';
    z_vector_KO = pinv(beta_prime_KO) * target_exp_valid_KO;

    bestThrIdx = find(endsWith(simulatedThresholds, best_threshold));
    thr = thresholds(bestThrIdx);

    z_bin_KO = double(z_vector_KO >= thr);
    y_tgt_KO = beta_matrix_KO * (TF_exp_KO.*z_bin_KO);
    y_rxn_KO = gene_to_rxn_map * y_tgt_KO;
    y_rxn_KO(isinf(y_rxn_KO)) = min(y_rxn_KO(~isinf(y_rxn_KO)));

    T = array2table(y_rxn_KO, ...
        'RowNames',modelWT_FBA.rxns,'VariableNames',{'WT'});

    for i = 1:length(z_bin_KO)
        zko = z_bin_KO;
        zko(i) = 0;
        yk  = gene_to_rxn_map*(beta_matrix_KO*(TF_exp_KO.*zko));
        yk(isinf(yk)) = min(yk(~isinf(yk)));
        T.(TF_list{i}) = yk;
    end

    Tnorm = T;
    for j_WT = 1:width(T)
        col = T{:, j_WT};
        tau_j = max(abs(col));
        if tau_j == 0
            Tnorm{:, j_WT} = 0;
        else
            Tnorm{:, j_WT} = col / tau_j;
        end
    end

    T = array2table(Tnorm{:,:}, ...
        'RowNames',modelWT_FBA.rxns,'VariableNames',T.Properties.VariableNames);

    y_KO_table = T;
    z_KO_table = z_bin_KO;

    KO_fluxes = struct();
    active_TFs_KO = TF_list(z_KO_table == 1);

    TF_growth_KO = table('Size', [length(active_TFs_KO)+1 1], ...
        'VariableTypes', {'double'}, ...
        'RowNames', ['WT'; active_TFs_KO], ...
        'VariableNames', {'Biomass'});
    active_TFs_KO = TF_growth_KO.Properties.RowNames;

    % --- Parfor Setup ---
    num_KOs = length(active_TFs_KO);
    biomass_par = zeros(num_KOs, 1);
    fluxes_par = cell(num_KOs, 1);

    y_KO_array = y_KO_table{:, :};
    rxnNames_y_KO = y_KO_table.Properties.RowNames;
    TF_names_y_KO = y_KO_table.Properties.VariableNames;

    % Restrict threads to prevent oversubscription across workers
    gurobi_params = struct();
    gurobi_params.OutputFlag = 0;
    gurobi_params.Threads = 1;

    % --- Parallel KO loop ---
    parfor i = 1:num_KOs
        TF = active_TFs_KO{i};
        modelKO = modelWT_FBA;

        if ismember(TF, modelKO.genes)
            modelKO = deleteModelGenes(modelKO, TF);
            [~, geneRxnIdx] = find(modelKO.rxnGeneMat(:, strcmp(modelKO.genes, TF)));
            for rxnIdx = geneRxnIdx'
                if ~isempty(modelKO.grRules{rxnIdx}) && contains(modelKO.grRules{rxnIdx}, TF)
                    modelKO.lb(rxnIdx) = 0;
                    modelKO.ub(rxnIdx) = 0;
                end
            end
        end

        S_KO = sparse(modelKO.S);
        [m_KO,n_KO] = size(S_KO);
        lb_base_KO = modelKO.lb;
        ub_base_KO = modelKO.ub;

        tf_idx_in_table = find(strcmp(TF_names_y_KO, TF), 1);
        y_rxn_KO_loop = y_KO_array(:, tf_idx_in_table);

        valid_KO = ~isnan(y_rxn_KO_loop);
        rxnIDs_KO = rxnNames_y_KO(valid_KO);
        y_rxn_KO_loop = y_rxn_KO_loop(valid_KO);
        y_rxn_KO_loop(~isfinite(y_rxn_KO_loop)) = 1;

        lbp_KO = lb_base_KO;
        ubp_KO = ub_base_KO;

        for j_KO = 1:numel(rxnIDs_KO)
            rxn_KO = rxnIDs_KO{j_KO};
            idx_KO = find(strcmp(modelKO.rxns, rxn_KO), 1);
            if isempty(idx_KO), continue; end

            if gene_to_rxn_map(idx_KO,:) == 0
                continue;
            end

            y_KO = y_rxn_KO_loop(j_KO);
            if ~isfinite(y_KO), y_KO = 1; end

            if abs(lb_base_KO(idx_KO)) > 1e-6 && abs(ub_base_KO(idx_KO)) > 1e-6
                lbp_KO(idx_KO) = max(lb_base_KO(idx_KO) * y_KO, -1000);
                ubp_KO(idx_KO) = min(ub_base_KO(idx_KO) * y_KO, 1000);
            elseif lb_base_KO(idx_KO) < -1e-6
                lbp_KO(idx_KO) = max(lb_base_KO(idx_KO) * y_KO, -1000);
            elseif ub_base_KO(idx_KO) > 1e-6
                ubp_KO(idx_KO) = min(ub_base_KO(idx_KO) * y_KO, 1000);
            end
        end

        lbp_KO(~isfinite(lbp_KO)) = lb_base_KO(~isfinite(lbp_KO));
        ubp_KO(~isfinite(ubp_KO)) = ub_base_KO(~isfinite(ubp_KO));

        infeasible_idx = find(lbp_KO > ubp_KO);
        if ~isempty(infeasible_idx)
            for k_inf = infeasible_idx'
                lbp_KO(k_inf) = min(lbp_KO(k_inf), ubp_KO(k_inf) - 1e-6);
                ubp_KO(k_inf) = max(lbp_KO(k_inf), ubp_KO(k_inf) + 1e-6);
            end
        end

        gp_KO = struct();
        gp_KO.A = sparse([
            S_KO,            sparse(m_KO,n_KO),    sparse(m_KO,n_KO);
           -speye(n_KO),    -speye(n_KO),          sparse(n_KO,n_KO);
            speye(n_KO),     sparse(n_KO,n_KO),   -speye(n_KO)
        ]);
        gp_KO.rhs        = [ zeros(m_KO,1); -lbp_KO; ubp_KO ];
        gp_KO.sense      = [ repmat('=', m_KO, 1); repmat('<', n_KO, 1); repmat('<', n_KO, 1) ];  % Use single chars: '>', '<'
        gp_KO.modelsense = 'min';
    
        gp_KO.obj = [ zeros(n_KO,1); WT_kappa*ones(n_KO,1); WT_kappa*ones(n_KO,1) ];
        gp_KO.obj(biomass_idx) = -mu;
    
        gp_KO.lb = [lb_base_KO; zeros(2*n_KO,1)];
        gp_KO.ub = [ub_base_KO; inf(2*n_KO,1)];
        gp_KO.varnames = [strcat('v_', modelKO.rxns); strcat('α_', modelKO.rxns); strcat('β_', modelKO.rxns)];

        result_KO = gurobi(gp_KO, gurobi_params);

        if isfield(result_KO, 'x')
            v_sol = result_KO.x(1:n_KO);
            biomass_par(i) = v_sol(biomass_idx);
            fluxes_par{i} = v_sol;
        else
            biomass_par(i) = 0;
            fluxes_par{i} = [];
        end
    end

    % Reconstruct the original data structures after parfor
    for i = 1:num_KOs
        TF = active_TFs_KO{i};
        TF_growth_KO{TF, 'Biomass'} = biomass_par(i);
        KO_fluxes.(TF) = fluxes_par{i};
    end

    % Evaluate predictions
    cutoff_range = 0.05:0.05:1.0;

    for cutoffVal = cutoff_range
        [lia, loc] = ismember(essentialityTable.Gene_ID, TF_growth_KO.Properties.RowNames);
        y_true     = essentialityTable.Essentiality_0_1(lia);
        if all(y_true == 0) || all(y_true == 1)
            continue;
        end

        biomassWT = TF_growth_KO.Biomass(1,1);
        cutoff = cutoffVal * biomassWT;

        pred_viable = TF_growth_KO >= cutoff;
        pred_viable = pred_viable{:,:};
        pred_viable = pred_viable(loc(lia));
        y_pred_ess  = double(~pred_viable);

        C = confusionmat(y_true, y_pred_ess);
        TN = C(1,1); FP = C(1,2);
        FN = C(2,1); TP = C(2,2);

        accuracy    = (TP+TN) / sum(C(:));
        sensitivity = TP / (TP + FN);
        specificity = TN / (TN + FP);

        scores = 1 - TF_growth_KO.Biomass(loc(lia))/biomassWT_FBA;
        [fpRate, tpRate, thr, AUC] = perfcurve(y_true, scores, 1);

        all_results.(selectedCondition).KO_growth  = TF_growth_KO.Biomass;
        all_results.(selectedCondition).KO_fluxes  = KO_fluxes;

        metrics_loop = table( ...
            string(selectedCondition), cutoffVal, accuracy, sensitivity, specificity, AUC, ...
            'VariableNames', {'Condition','cutoffVal','Accuracy','Sensitivity','Specificity','AUC'} ...
        );
        metrics_loop.tpRate{1} = tpRate;
        metrics_loop.fpRate{1} = fpRate;
        metrics_loop.thr{1}    = thr;
        metrics_loop.C{1}      = C;

        metrics_all = [metrics_all; metrics_loop];
    end
end

metrics_all = rmmissing(metrics_all);

%% Final global metrics across all conditions
global_acc  = mean(metrics_all.Accuracy);
global_sens = mean(metrics_all.Sensitivity);
global_spec = mean(metrics_all.Specificity);

fprintf('\n[%s] === Global metrics across all conditions ===', ...
    datetime('now','Format','HH:mm:ss'))

fprintf('\n[%s] Accuracy = %.3f | Sensitivity = %.3f | Specificity = %.3f\n', ...
    datetime('now','Format','HH:mm:ss'), ...
    global_acc, global_sens, global_spec);

metrics_OCELOT_yeast = metrics_all;
TF_growth_KO_OCELOT_yeast = TF_growth_KO;

save('Sce_essentiality.mat', 'metrics_OCELOT_yeast', 'TF_growth_KO_OCELOT_yeast', '-v7.3');
fprintf('\n[%s] === FINISHED (results saved to disk) ===\n', ...
    datetime('now','Format','HH:mm:ss'))

%%
toc
