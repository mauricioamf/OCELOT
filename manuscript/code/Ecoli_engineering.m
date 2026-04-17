%% Initial configuration
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);
clc;clear

tic

%% Load metabolic model and other data

% Load metabolic model
load('../models/iML1515.mat');

% Load expression data
expression_data = readtable('../data/Ecoli_count_matrix.csv', 'ReadRowNames', true, 'VariableNamingRule', 'preserve');

% Load beta matrix as dataframe
beta_df = readtable('../data/Ecoli_GRN_lasso.csv', 'ReadRowNames', true);

%% Set parameters
% General parameters
warning("on");
showPlots = 0; % 0 = no plots, 1 = confusion matrix and ROC only, 2 = all

% WT optimization steps
selectedReference = "control__wt_glc"; % Cluster_1 or control__wt_glc or MG1655
biomass_rxn = 'BIOMASS_Ec_iML1515_core_75p37M';
glucoseUptake = -10;
mu = 1;                       % only applicable if runWToptimization = true
slack_max = 1000;
saveWToptimization = false;
thresholdSelection = 'Composite';  % 'RMSE' or 'Growth'
max_infeas   = 3;  % stop after n consecutive infeasible thresholds

% Engineering steps
prod_rxn = 'EX_ac_e'; % EX_ac_e EX_etoh_e EX_for_e EX_lac__D_e EX_succ_e  ONE AT A TIME
max_iterations = 12;
max_perturbations = 5;
nKOsLimit = 20;

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
    % avg_expr_data{:,c} = min(expression_data{:,idx}, [], 2);  
end

%% Optimizations for WT
% save('checkpoint_WT.mat')
% clc;clear
% load('checkpoint_WT.mat')

% FBA for WT growth rate
modelWT_FBA = model;

exchangeIndex = find(contains(modelWT_FBA.rxnNames, "exchange"));
exchangeIDs = modelWT_FBA.rxns(exchangeIndex);
exchangeIDs(end) = [];
biomass_idx = find(strcmp(modelWT_FBA.rxns, biomass_rxn));

M9_components = ["EX_pi_e", "EX_co2_e", "EX_fe3_e", "EX_h_e", ...
    "EX_mn2_e", "EX_fe2_e", "EX_zn2_e", "EX_mg2_e", ...
    "EX_ca2_e", "EX_ni2_e", "EX_cu2_e", "EX_sel_e", ...
    "EX_cobalt2_e", "EX_h2o_e", "EX_mobd_e", "EX_so4_e", ...
    "EX_nh4_e", "EX_k_e", "EX_na1_e", "EX_cl_e", ...
    "EX_o2_e", "EX_tungs_e", "EX_slnt_e"];

modelWT_FBA = changeRxnBounds(modelWT_FBA, exchangeIDs, 0, 'l');
modelWT_FBA = changeRxnBounds(modelWT_FBA, M9_components, -1000, 'l');
modelWT_FBA = changeObjective(modelWT_FBA, modelWT_FBA.rxns{biomass_idx});

modelWT_FBA = changeRxnBounds(modelWT_FBA, 'EX_glc__D_e', glucoseUptake, 'b'); %D-glucose
FBAsol = optimizeCbModel(modelWT_FBA);
% printFluxes(modelWT_FBA, FBAsol.x, true);
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
tf_vals    = TF_exp;
thresholds = unique(sort(tf_vals));   % one threshold per TF expression level

% Output structs
y_all_struct   = struct();  % holds y_all_tau per threshold
z_struct       = struct();  % holds z_vec per threshold
z_struct_cont  = struct();  % holds z_vec per threshold

% Loop through thresholds
for t = 1:numel(thresholds)
    thr     = thresholds(t);
    z_bin = double(TF_exp >= thr);
    z_cont = TF_exp;
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
    z_cont_struct.(name) = y_tgt;
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
    fprintf('\n--- WT PROM for threshold: %s ---\n', threshold_name);

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
        target_exp_hat = z_cont_struct.(threshold_name);
    
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
                % For reversible reactions
                lbp_WT(idx_WT) = max(lb_base_WT(idx_WT) * y_WT, -1000);
                ubp_WT(idx_WT) = min(ub_base_WT(idx_WT) * y_WT, 1000);
            elseif lb_base_WT(idx_WT) < -1e-6
                % For uptake reactions (negative lb)
                lbp_WT(idx_WT) = max(lb_base_WT(idx_WT) * y_WT, -1000);
            elseif ub_base_WT(idx_WT) > 1e-6  
                % For production reactions (positive ub)
                ubp_WT(idx_WT) = min(ub_base_WT(idx_WT) * y_WT, 1000);
            end

        end

        lbp_WT(~isfinite(lbp_WT)) = lb_base_WT(~isfinite(lbp_WT));
        ubp_WT(~isfinite(ubp_WT)) = ub_base_WT(~isfinite(ubp_WT));

        infeasible_idx = find(lbp_WT > ubp_WT);
        if ~isempty(infeasible_idx)
            fprintf('Found %d reactions with lb > ub after scaling — fixing them.\n', numel(infeasible_idx));
            for idx_WT = infeasible_idx'
                lbp_WT(idx_WT) = min(lbp_WT(idx_WT), ubp_WT(idx_WT) - 1e-6);  % or set lb = ub = 0
                ubp_WT(idx_WT) = max(lbp_WT(idx_WT), ubp_WT(idx_WT) + 1e-6);
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

                    fprintf('→ predicted biomass = %.4f, slacks = %d, kappa = %.3e', ...
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
                fprintf('%s → No solution at kappa = %.3e\n', threshold_name, kappa_WT);
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
            fprintf('\nStopping WT simulations after %d consecutive infeasible thresholds.\n', max_infeas);
            % break;
        end
    end
end

%% Compute RMSE for calculated target_exp_hat
simulatedThresholds = fieldnames(WT_results);

rmse_table = table('Size', [length(simulatedThresholds), 4], ...
    'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Threshold', 'RMSE', 'R2', 'Growth'});

% Loop through thresholds
for th = 1:length(simulatedThresholds)
    threshold_name = simulatedThresholds{th};
    predicted_tbl = WT_results.(threshold_name).target_expr;
    predicted_vec = predicted_tbl{:,1};  % predicted target expression (G×1)
    target_exp_valid = target_exp(validTargets);

    predicted_vec = log(abs(predicted_vec));
    target_exp_valid = log(abs(target_exp_valid));

    % Calculate RMSE
    residual = target_exp_valid - predicted_vec;

    goodIdx = isfinite(residual);
    target_exp_valid = target_exp_valid(goodIdx);
    predicted_vec = predicted_vec(goodIdx);

    residual = target_exp_valid - predicted_vec;
    rmse_value = sqrt(mean(residual.^2));

    % Calculate R²
    ss_res = sum((predicted_vec - target_exp_valid).^2);
    ss_tot = sum((target_exp_valid - mean(target_exp_valid)).^2);
    r2 = 1 - (ss_res / ss_tot);

    % Store
    rmse_table.Threshold(th) = threshold_name;
    rmse_table.RMSE(th) = rmse_value;
    rmse_table.R2(th) = r2;
    rmse_table.Growth(th) = WT_results.(threshold_name).biomass;
end

rmse_table.Discriminability = zeros(height(rmse_table), 1);

for th = 1:length(simulatedThresholds)
    threshold_name = simulatedThresholds{th};

    % Recover z_bin for this threshold
    z_bin_th = z_struct.(threshold_name);
    n_active = sum(z_bin_th);

    % If all or none are active, discriminability is zero — skip
    if n_active == 0 || n_active == numel(z_bin_th)
        rmse_table.Discriminability(th) = 0;
        continue;
    end

    % Baseline reaction scaling with all active TFs
    y_base = gene_to_rxn_map * (beta_matrix * (TF_exp .* z_bin_th));

    % For each TF, compute how much its KO shifts the reaction scaling
    delta_norms = zeros(numel(TF_list), 1);
    for i = 1:numel(TF_list)
        if z_bin_th(i) == 0
            % Already inactive — KO has no additional effect
            delta_norms(i) = 0;
            continue;
        end
        zko = z_bin_th;
        zko(i) = 0;
        y_ko = gene_to_rxn_map * (beta_matrix * (TF_exp .* zko));
        delta_norms(i) = norm(y_base - y_ko);
    end

    % Discriminability = coefficient of variation of KO impacts
    % (std/mean): high when TFs have heterogeneous effects on flux
    active_deltas = delta_norms(z_bin_th == 1);
    if mean(active_deltas) > 0
        rmse_table.Discriminability(th) = std(active_deltas) / mean(active_deltas);
    else
        rmse_table.Discriminability(th) = 0;
    end
end

% Select best threshold
valid_idx = ~isnan(rmse_table.RMSE) & ~isnan(rmse_table.Growth) & rmse_table.Growth > 0;

switch thresholdSelection
    case 'Growth'
        [~, best_local_idx] = max(rmse_table.Growth(valid_idx));

    case 'RMSE'
        [~, best_local_idx] = min(rmse_table.RMSE(valid_idx));

    case 'Discriminability'
        [~, best_local_idx] = max(rmse_table.Discriminability(valid_idx));

    case 'Composite'
        rmse_valid   = rmse_table.RMSE(valid_idx);
        disc_valid   = rmse_table.Discriminability(valid_idx);
        growth_valid = rmse_table.Growth(valid_idx);
        r2_valid     = rmse_table.R2(valid_idx);
    
        % Normalise to [0,1]; RMSE inverted (lower is better), R2 clipped to [0,1]
        rmse_norm   = 1 - (rmse_valid - min(rmse_valid))     / (max(rmse_valid)   - min(rmse_valid)   + eps);
        disc_norm   =     (disc_valid - min(disc_valid))     / (max(disc_valid)   - min(disc_valid)   + eps);
        growth_norm =     (growth_valid - min(growth_valid)) / (max(growth_valid) - min(growth_valid) + eps);
        r2_norm     = max(r2_valid, 0);  % R2 can be negative (worse than mean); clip to 0
    
        composite = 0.4 * disc_norm + 0.3 * r2_norm + 0.2 * rmse_norm + 0.1 * growth_norm;
        [~, best_local_idx] = max(composite);
end

best_idx_all = find(valid_idx);
best_idx     = best_idx_all(best_local_idx);
best_threshold = rmse_table.Threshold(best_idx);

fprintf('\nSelected threshold: %s\n', best_threshold);
fprintf('  RMSE=%.4f | R2=%.4f | Growth=%.4f | Discriminability=%.4f\n', ...
    rmse_table.RMSE(best_idx), ...
    rmse_table.R2(best_idx), ...
    rmse_table.Growth(best_idx), ...
    rmse_table.Discriminability(best_idx));

%% Engineering pipeline steps 1 and 2

% Step 1: maximize biomass
results_best_WT = WT_results.(best_threshold);

% Step 2: maximize product
prod_idx = find(strcmp(model.rxns,prod_rxn));
prod_name = model.rxnNames(prod_idx);
% flex = 0.05;

modelProd_FBA = modelWT_FBA;
modelProd_FBA.c(:) = 0;
modelProd_FBA.c(prod_idx) = 1;
modelProd_FBA.lb(prod_idx) = -1000;
modelProd_FBA.ub(prod_idx) = 1000;
modelProd_FBA.lb(biomass_idx) = biomassWT_FBA * 0.5;
solProd_FBA = optimizeCbModel(modelProd_FBA);
% printFluxes(modelProd_FBA, solProd_FBA.x, true);

%% Engineering pipeline step 3: minimize z rearrangement with integer cuts

% save('checkpoint_step3.mat')
% clear
% load('checkpoint_step3.mat')

% Compute WT reaction activity baseline for ratio-based bound scaling
z_best    = z_struct.(best_threshold);
y_WT_base = gene_to_rxn_map * (beta_matrix * (TF_exp .* z_best));
y_WT_base(isinf(y_WT_base)) = min(y_WT_base(~isinf(y_WT_base)));

results_all = {};

z0 = z_struct.(best_threshold);

[nMets, nRxns] = size(modelWT_FBA.S);
nTF = numel(TF_list);

% Build regulatory influence matrix
upsilon = gene_to_rxn_map * beta_prime;   % nRxns × nTF
upsilon = upsilon ./ max(upsilon,[],'all');

% Offsets in x = [ v(1..nRxns); α(1..nRxns); β(1..nRxns); z(1..nTF) ]
v0     = 0;
a0     = v0 + nRxns;
b0     = a0  + nRxns;
z0_off = b0  + nRxns;
nVars  = z0_off + nTF;

% Variable types
vtype = [ repmat('C', nRxns,1);
          repmat('C', nRxns,1);
          repmat('C', nRxns,1);
          repmat('B', nTF,1) ];

% Gurobi settings
params = struct();
params.OutputFlag = 0;
params.TimeLimit = 600;
params.MIPGap = 0.01;
params.MIPFocus = 1;

flex_range = 0.5;
kappa_range = WT_results.(best_threshold).kappa;  

% Initialize a cell array to store integer cuts
integer_cuts = {};  % Each cut will be a struct with 'z_pattern', 'sum_value'

for flex = flex_range

    fb = flex;  % Same flex for biomass
    fp = flex;  % Same flex for product

    for kappa = kappa_range

        min_biomass_step3 = biomassWT_FBA * fb;
        min_product_step3 = solProd_FBA.f * fp;

        % Base bounds
        lb_base = modelWT_FBA.lb;
        ub_base = modelWT_FBA.ub;
        lb_base(biomass_idx) = min_biomass_step3;
        lb_base(prod_idx)    = min_product_step3;

        % Build constraint arrays
        A   = sparse(0, nVars);
        rhs = [];
        sen = '';
        
        % Steady‐state: S*v = 0
        A   = [A; sparse(modelWT_FBA.S), sparse(nMets, nVars - nRxns)];
        rhs = [rhs; zeros(nMets,1)];
        sen = [sen; repmat('=', nMets, 1)];

        % PROM slacks for each regulated reaction
        for i = 1:nRxns
            if any(upsilon(i,:))
                % (a) v_i + α_i ≥ lb_i * (D(i,:)·z)
                row = zeros(1, nVars);
                row(v0 + i) =  1;
                row(a0 + i) =  1;
                row(z0_off+1 : z0_off+nTF) = -lb_base(i) * upsilon(i,:);
                A   = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '>'];
        
                % (b) v_i - β_i ≤ ub_i * (D(i,:)·z)
                row = zeros(1, nVars);
                row(v0 + i) =  1;
                row(b0 + i) = -1;
                row(z0_off+1 : z0_off+nTF) = ub_base(i) * upsilon(i,:);
                A   = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '<'];
            end
        end

        perturb_row = zeros(1, nVars);
        for j = 1:nTF
            if z0(j) == 1
                % For wild-type active TFs: (1 - z_j) contributes to KO count
                perturb_row(z0_off + j) = -1;  % -z_j
            else
                % For wild-type inactive TFs: z_j contributes to activation count
                perturb_row(z0_off + j) = 1;   % z_j
            end
        end

        % The RHS is max_perturbations - sum(z0==1)
        % Because sum_{j: z0(j)=1} (1 - z_j) = sum(z0==1) - sum_{j: z0(j)=1} z_j
        rhs_perturb = max_perturbations - sum(z0 == 1);

        A = [A; perturb_row];
        rhs = [rhs; rhs_perturb];
        sen = [sen; '<'];
       
        % Variable bounds and types
        lb = [ lb_base;         % v
               zeros(2*nRxns,1);% α,β
               zeros(nTF,1) ];  % z
        
        ub = [ ub_base;         % v
               inf(2*nRxns,1);  % α,β
               ones(nTF,1) ];   % z
  
        % Reset integer cuts for this parameter combination
        current_integer_cuts = {};
        
        % Iterative integer cut strategy
        for iteration = 1:max_iterations
            
            % Add existing integer cuts to the model
            A_with_cuts = A;
            rhs_with_cuts = rhs;
            sen_with_cuts = sen;
            
            % Add all accumulated integer cuts
            for cut_idx = 1:length(integer_cuts)
                cut = integer_cuts{cut_idx};
                row = zeros(1, nVars);
                
                % Create cut: sum of (z_i where pattern_i=1) + sum of (1-z_i where pattern_i=0) ≤ nTF - 2
                % This excludes the exact pattern
                for j = 1:nTF
                    if cut.z_pattern(j) == 1
                        row(z0_off + j) = 1;
                    else
                        row(z0_off + j) = -1;
                    end
                end
                
                A_with_cuts = [A_with_cuts; row];
                rhs_with_cuts = [rhs_with_cuts; cut.sum_value - 2];  % Make it <= sum_value - 1
                sen_with_cuts = [sen_with_cuts; '<'];
            end
                       
            % Objective: reward preserving z₀, penalize slacks
            obj = zeros(nVars, 1);
            
            for j = 1:nTF
                if z0(j) == 1
                    obj(z0_off + j) =  +1;   % reward z=1
                else
                    obj(z0_off + j) =  -1;   % reward z=0  (maximize -z)
                end
            end
            
            obj(a0+1 : b0) = -kappa;       % penalize α
            obj(b0+1 : z0_off) = -kappa;   % penalize β
            
            model_step3 = struct();
            model_step3.A = A_with_cuts;
            model_step3.rhs = rhs_with_cuts;
            model_step3.sense = sen_with_cuts;
            model_step3.lb = lb;
            model_step3.ub = ub;
            model_step3.vtype = vtype;
            model_step3.obj = obj;
            model_step3.modelsense = 'max';
        
            result_step3 = gurobi(model_step3, params);
            
            if strcmp(result_step3.status, 'OPTIMAL')
                z_sol = result_step3.x(z0_off+1 : z0_off+nTF);
                num_ko = sum(z0 == 1 & z_sol < 0.5);
        
                v_sol      = result_step3.x(1:nRxns);
                alpha_vals = result_step3.x(a0+1 : b0);
                beta_vals  = result_step3.x(b0+1 : z0_off);
    
                knockout_idx = find(z0 == 1 & z_sol < 0.5);
                knocked_out_TFs = TF_list(knockout_idx);
                activated_idx = find(z0 == 0 & z_sol > 0.5);
                activated_TFs = TF_list(activated_idx);
    
                biomass_eng = v_sol(biomass_idx);
                product_eng = v_sol(prod_idx);
                num_slacks_used = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);
        
                % Collect into results
                results_all(end+1,:) = {fb, fp, biomass_eng, product_eng, ...
                                        knocked_out_TFs, activated_TFs, ...
                                        numel(knocked_out_TFs), numel(activated_TFs), ...
                                        num_slacks_used, kappa, z0, z_sol, iteration, 'Original'};
                
                % Store the z pattern to prevent it from appearing again
                cut_pattern = struct();
                cut_pattern.z_pattern = round(z_sol);  % Round to 0 or 1
                cut_pattern.sum_value = sum(cut_pattern.z_pattern);  % Sum of z_i where pattern_i=1
                
                % Check if this pattern is already in our cuts
                is_new_pattern = true;
                for cut_idx = 1:length(integer_cuts)
                    if isequal(integer_cuts{cut_idx}.z_pattern, cut_pattern.z_pattern)
                        is_new_pattern = false;
                        break;
                    end
                end
                
                if is_new_pattern
                    integer_cuts{end+1} = cut_pattern;
                    current_integer_cuts{end+1} = cut_pattern;
                end
                
            else
                % If no solution found, break the iteration loop
                if iteration == 1
                    results_all(end+1,:) = {fb, fp, NaN, NaN, {}, {}, 0, 0, NaN, kappa, z0, NaN, 0, 'Original'};
                end
                break;
            end
            
           
        end  % End iteration loop
    end  % End kappa loop
end  % End fb loop

% Save results
resultsTable = cell2table(results_all, ...
    'VariableNames', {'BiomassFlex','ProductFlex','Biomass','Product', ...
                      'Knockouts','Activated','NumKO','NumActivated','SlacksUsed','Kappa', 'z0', 'z_sol', 'Iteration', 'Strategy'});

resultsTable = resultsTable(:, {'Strategy','Iteration','BiomassFlex','ProductFlex','Kappa','Biomass','Product', ...
                               'Knockouts','Activated','NumKO','NumActivated','SlacksUsed','z0','z_sol'});

%% Checkpoint
% save('checkpoint_integetcuts.mat')
% clc;clear
% load('checkpoint_integetcuts.mat')

%% Process results
% rows with NaNs in any numeric column
numericData = table2array(resultsTable(:, vartype('numeric')));  % selects numeric variables only
rowHasNaN = any(isnan(numericData), 2);
% rows where both NumKO and NumActivated are zero (no perturbation)
rowNoPerturb = (resultsTable.NumKO == 0) & (resultsTable.NumActivated == 0);
% rows where KO or Activation affects ALL TFs (unrealistic)
rowAllKOorAllAct = (resultsTable.NumKO == nTF) | (resultsTable.NumActivated == nTF);
% rows where z_sol is all zeroes
rowZallZeroes = cellfun(@(v) all(v == 0), resultsTable.z_sol);
% combine masks and filter
badRows = rowHasNaN | rowNoPerturb | rowAllKOorAllAct | rowZallZeroes;
fprintf('\nFiltering %d / %d rows (NaN / no-perturb / all-TF)...\n', sum(badRows), height(resultsTable));
resultsTableFilter = resultsTable(~badRows, :);

%% Validation: predicting growth and product with predicted z_sol
fprintf('\n--- Validation: predicting growth and product with predicted z_sol---\n');

% Preallocate results
% KO_names = TF_list(knockout_idx);
nSols         = height(resultsTableFilter);
KO_names      = cell(nSols,1);
act_names     = cell(nSols,1);
strategy_val  = cell(nSols,1);
iteration_val = NaN(nSols,1);

growth_val_vbio = NaN(nSols,1);
prod_val_vbio   = NaN(nSols,1);
slack_val_vbio  = NaN(nSols,1);

growth_val_prod = NaN(nSols,1);
prod_val_prod   = NaN(nSols,1);
slack_val_prod  = NaN(nSols,1);

% Use the same D matrix as in engineering step
upsilon_val = gene_to_rxn_map * beta_prime;
upsilon_val = upsilon_val ./ max(upsilon_val(:));  % Same normalization as engineering

% Build the problem with regulatory constraints (same as engineering)
[nMets_val, nRxns_val] = size(modelWT_FBA.S);
nTF_val = numel(TF_list);

% Variable indices
v0_val    = 0;
a0_val    = v0_val + nRxns_val;
b0_val    = a0_val + nRxns_val;
z0_val    = b0_val + nRxns_val;
nVars_val = z0_val + nTF_val;

% Variable types (z is fixed, not binary in validation)
vtype_val = [repmat('C', nRxns_val,1);    % v
             repmat('C', nRxns_val,1);    % α
             repmat('C', nRxns_val,1);    % β
             repmat('B', nTF_val,1)];     % z (continuous, fixed)

% Build constraint arrays
A_val = sparse(0, nVars_val);
rhs_val = [];
sen_val = '';

% Steady-state: S*v = 0
A_val = [A_val; sparse(modelWT_FBA.S), sparse(nMets_val, nVars_val - nRxns_val)];
rhs_val = [rhs_val; zeros(nMets_val,1)];
sen_val = [sen_val; repmat('=', nMets_val, 1)];

% Regulatory constraints
for i = 1:nRxns_val
    if any(upsilon_val(i,:))
        % alpha
        row = zeros(1, nVars_val);
        row(v0_val + i) =  1;
        row(a0_val + i) =  1;
        row(z0_val+1 : z0_val+nTF_val) = -lb_base(i) * upsilon_val(i,:);
        A_val = [A_val; row];
        rhs_val = [rhs_val; 0];
        sen_val = [sen_val; '>'];
        
        % beta
        row = zeros(1, nVars_val);
        row(v0_val + i) =  1;
        row(b0_val + i) = -1;
        row(z0_val+1 : z0_val+nTF_val) = ub_base(i) * upsilon_val(i,:);
        A_val = [A_val; row];
        rhs_val = [rhs_val; 0];
        sen_val = [sen_val; '<'];
    end
end

% Use the engineered z vector
for z = 1:height(resultsTableFilter)

    KO_names{z}      = resultsTableFilter.Knockouts{z,1};
    act_names{z}    = resultsTableFilter.Activated{z,1};
    strategy_val{z}  = resultsTableFilter.Strategy{z,1};
    iteration_val(z) = resultsTableFilter.Iteration(z,1);
   
    z_val = resultsTableFilter.z_sol{z,1};
    kappa_val = resultsTableFilter.Kappa(z,1);

    % if nSols < 10
        fprintf('Validating engineered solution %d of %d\n', z, nSols);
    % end
    
    % Variable bounds
    lb_val = [modelWT_FBA.lb;                  % v
              zeros(2*nRxns_val,1);            % α, β 
              z_val];                          % z
    
    ub_val = [modelWT_FBA.ub;                  % v
              inf(2*nRxns_val,1);              % α, β
              z_val];                          % z
    
    % Add biomass/product flexibility constraints
    lb_val(v0_val + biomass_idx) = resultsTableFilter.Biomass(z);
    lb_val(v0_val + prod_idx)    = resultsTableFilter.Product(z);
    
    % Optimize for biomass
    model_val_bio = struct();
    model_val_bio.A = A_val;
    model_val_bio.rhs = rhs_val;
    model_val_bio.sense = sen_val;
    model_val_bio.lb = lb_val;
    model_val_bio.ub = ub_val;
    model_val_bio.vtype = vtype_val;
    model_val_bio.modelsense = 'max';
    
    model_val_bio.obj = zeros(nVars_val, 1);
    model_val_bio.obj(v0_val + biomass_idx) = mu;
    % model_val_bio.obj(a0_val+1 : b0_val) = -kappa_val;
    % model_val_bio.obj(b0_val+1 : z0_val) = -kappa_val;
    
    params.OutputFlag = 0;
    result_val_bio = gurobi(model_val_bio, params);
    
    if strcmp(result_val_bio.status, 'OPTIMAL')
        v_sol = result_val_bio.x(1:nRxns_val);
        alpha_vals = result_val_bio.x(nRxns_val+1:2*nRxns_val);
        beta_vals = result_val_bio.x(2*nRxns_val+1:3*nRxns_val);
        
        biomass_val = v_sol(biomass_idx);
        product_val = v_sol(prod_idx);
        num_slacks = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

        growth_val_vbio(z) = biomass_val;
        prod_val_vbio(z) = product_val;
        slack_val_vbio(z) = num_slacks;
        
        fprintf('  → Biomass = %.4f (%.2f%% WT), %s = %.4f (%.2f%% of maximum), slacks = %d\n', ...
                    biomass_val, 100*biomass_val/biomassWT_FBA, ...
                    prod_name{1,1}, product_val, 100*product_val/solProd_FBA.f, ...
                    num_slacks);
    
    end
    
    % Solve for product
    model_val_prod = model_val_bio;
    model_val_prod.obj = zeros(nVars_val, 1);
    model_val_prod.obj(v0_val + prod_idx) = 1;
    % model_val_prod.obj(a0_val+1 : b0_val) = -kappa_val;
    % model_val_prod.obj(b0_val+1 : z0_val) = -kappa_val;
   
    result_val_prod = gurobi(model_val_prod, params);
    
    if strcmp(result_val_prod.status, 'OPTIMAL')
        v_sol = result_val_prod.x(1:nRxns_val);
        biomass_val = v_sol(biomass_idx);
        product_val = v_sol(prod_idx);

        alpha_vals = result_val_prod.x(nRxns_val+1:2*nRxns_val);
        beta_vals = result_val_prod.x(2*nRxns_val+1:3*nRxns_val);
        num_slacks = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

        growth_val_prod(z) = biomass_val;
        prod_val_prod(z) = product_val;
        slack_val_prod(z) = num_slacks;
        
        fprintf('  → Biomass = %.4f (%.2f%% WT), %s = %.4f (%.2f%% of maximum), slacks = %d\n', ...
                    biomass_val, 100*biomass_val/biomassWT_FBA, ...
                    prod_name{1,1}, product_val, 100*product_val/solProd_FBA.f, ...
                    num_slacks);
    end

end

%% Process output for exporting
final_results = table('Size', [height(resultsTableFilter), 10], ...
    'VariableTypes', {'cell','cell','cell','cell', ...
                      'double','double','double','double' ...
                      'double','double'}, ...
    'VariableNames', {'Strategy','Iteration','Knockouts','Activations', ...
                      'Max_growth','Min_product','Min_growth','Max_product', ...
                      'Slack_growth','Slack_product'});

final_results.Strategy      = strategy_val;
final_results.Iteration     = iteration_val;
final_results.Knockouts     = KO_names;
final_results.Activations   = act_names;
final_results.Max_growth    = growth_val_vbio;
final_results.Min_product   = prod_val_vbio;
final_results.Min_growth    = growth_val_prod;
final_results.Max_product   = prod_val_prod;
final_results.Slack_growth  = slack_val_vbio;
final_results.Slack_product = slack_val_prod;

% Format gene names
KOconcatName  = strings(height(final_results),1);
ACTconcatName = strings(height(final_results),1);

for i = 1:height(final_results)
    % Knockouts
    ko = final_results.Knockouts{i}; 
    if isempty(ko)
        KOconcatName(i) = "";
    else
        KOconcatName(i) = strjoin(string(ko), ", ");
    end
    % Activations
    act = final_results.Activations{i};
    if isempty(act)
        ACTconcatName(i) = "";
    else
        ACTconcatName(i) = strjoin(string(act), ", ");
    end
end
final_results.Knockouts   = KOconcatName;
final_results.Activations = ACTconcatName;

% Retrieve list of genes
nRows = height(final_results);
rxnNames = modelWT_FBA.rxns; 
AffectedRxnNames = strings(nRows,1);
final_results.AffectedReactions = strings(height(final_results),1);

for i = 1:nRows

    koTFs  = KO_names{i};
    actTFs = act_names{i};

    TFs_row = {};
    if ~isempty(koTFs)
        TFs_row = [TFs_row, koTFs];
    end
    if ~isempty(actTFs)
        TFs_row = [TFs_row, actTFs];
    end

    % No TFs → no affected reactions
    if isempty(TFs_row)
        AffectedRxnNames(i) = "";
        continue
    end

    % Map TF names to indices
    tf_idx = find(ismember(TF_list, TFs_row));

    if isempty(tf_idx)
        AffectedRxnNames(i) = "";
        continue
    end

    % Find affected reactions
    affected_rxn_idx = find(any(upsilon_val(:, tf_idx) ~= 0, 2));

    if isempty(affected_rxn_idx)
        AffectedRxnNames(i) = "";
    else
        AffectedRxnNames(i) = strjoin(string(rxnNames(affected_rxn_idx)), ", ");
    end
end

% Store in table
final_results.AffectedReactions = AffectedRxnNames;

% Add WT results
WTmaxBiomass = biomassWT_FBA;
WTminProduct = FBAsol.x(prod_idx);
WTminBiomass = solProd_FBA.x(biomass_idx);
WTmaxProduct = solProd_FBA.x(prod_idx);

WTresults = {'WT', 0, '', '', WTmaxBiomass, WTminProduct, WTminBiomass, WTmaxProduct,  0, 0, ''};

final_results = [WTresults; final_results];

% Export
metName = model.metNames(find(model.S(:,prod_idx)));
exportFilename = strcat('Ecoli_strategies_', metName, '.xlsx');
% writetable(final_results, exportFilename{1,1}, "FileType", "spreadsheet")

%% Knockout each predicted TF individually (single KO)
fprintf('\n--- Validating Single Knockouts ---\n');

if ~isstring(final_results.Knockouts)
    knockout_strings = string(final_results.Knockouts);
else
    knockout_strings = final_results.Knockouts;
end

% Split each string by ", " and flatten the results to get all tested TFs
all_genes = [];
for i = 1:length(knockout_strings)
    if knockout_strings(i) ~= ""
        split_genes = strsplit(knockout_strings(i), ', ');
        split_genes = strtrim(split_genes);
        all_genes = [all_genes; split_genes(:)];
    end
end

% Get unique genes
unique_genes = unique(all_genes);
knockout_idx = find(contains(TF_list, unique_genes));

% Preallocate results
KO_names = unique_genes;
nKOs = numel(KO_names);
growth_sKO_vbio  = NaN(nKOs,1);
prod_sKO_vbio    = NaN(nKOs,1);
slack_sKO_vbio   = NaN(nKOs,1);
growth_sKO_prod  = NaN(nKOs,1);
prod_sKO_prod    = NaN(nKOs,1);
slack_sKO_prod   = NaN(nKOs,1);

if nKOs >= nKOsLimit
    fprintf('Too many single KOs (nKOs = %d), suppressing output\n', nKOs);
end

% Retrieve the optimal kappa and z from the WT step
kappa_val = WT_results.(best_threshold).kappa;
z_base = z_struct.(best_threshold);

[nMets_sKO, nRxns_sKO] = size(modelWT_FBA.S);
v0_sKO = 0;
a0_sKO = nRxns_sKO;
b0_sKO = 2*nRxns_sKO;
nVars_sKO = 3 * nRxns_sKO; 

vtype_sKO = repmat('C', nVars_sKO, 1);

for k = 1:nKOs
    tf_idx = knockout_idx(k);
    if nKOs < nKOsLimit
        fprintf('Validating KO of TF %s (%d of %d)\n', KO_names{k}, k, nKOs);
    end

    % Build z vector with only this TF knocked out
    z_sKO = z_base;
    z_sKO(tf_idx) = 0;

    % Compute normalized y
    TF_activity_sKO = TF_exp .* z_sKO;  
    y_rxn_sKO = beta_matrix * TF_activity_sKO;
    y_rxn_sKO = gene_to_rxn_map * y_rxn_sKO;
    y_rxn_sKO = y_rxn_sKO ./ max(y_rxn_sKO,[],'all');

    % Handle numerical issues
    y_rxn_sKO(isnan(y_rxn_sKO)) = 1;
    y_rxn_sKO(~isfinite(y_rxn_sKO)) = 1;

    % Build engineering formulation with FIXED z and kappa penalty
    A_sKO = sparse(0, nVars_sKO);
    rhs_sKO = [];
    sen_sKO = '';

    % Steady-state: S*v = 0
    A_sKO   = [A_sKO; sparse(modelWT_FBA.S), sparse(nMets_sKO, 2*nRxns_sKO)];
    rhs_sKO = [rhs_sKO; zeros(nMets_sKO,1)];
    sen_sKO = [sen_sKO; repmat('=', nMets_sKO, 1)];

    % Regulatory bounds via slacks
    for i = 1:nRxns_sKO
        y_i = y_rxn_sKO(i);

        % Apply constraint only if reaction is regulated and scaled
        if any(gene_to_rxn_map(i,:)) && y_i < 0.99 
            % Lower bound constraint: v_i + α_i ≥ lb_i * y_i
            row_lb = zeros(1, nVars_sKO);
            row_lb(v0_sKO + i) = 1;
            row_lb(a0_sKO + i) = 1;
            A_sKO = [A_sKO; row_lb];
            rhs_sKO = [rhs_sKO; modelWT_FBA.lb(i) * y_i];
            sen_sKO = [sen_sKO; '>'];

            % Upper bound constraint: v_i - β_i ≤ ub_i * y_i
            row_ub = zeros(1, nVars_sKO);
            row_ub(v0_sKO + i) = 1;
            row_ub(b0_sKO + i) = -1;
            A_sKO = [A_sKO; row_ub];
            rhs_sKO = [rhs_sKO; modelWT_FBA.ub(i) * y_i];
            sen_sKO = [sen_sKO; '<'];
        end
    end  

    % Variable bounds (Absolute limits)
    lb_sKO = [modelWT_FBA.lb; zeros(2*nRxns_sKO,1)]; % v bounds, non-negative slacks
    ub_sKO = [modelWT_FBA.ub; inf(2*nRxns_sKO,1)];   % v bounds, inf slacks

    % Apply minimum biomass and product requirements (if applicable)
    lb_sKO(v0_sKO + biomass_idx) = biomassWT_FBA * 0.5;
    ub_sKO(v0_sKO + biomass_idx) = 1000;

    lb_sKO(v0_sKO + prod_idx) = solProd_FBA.f * 0.5;
    ub_sKO(v0_sKO + prod_idx) = 1000;

    % Build base Gurobi model structure
    model_sKO = struct();
    model_sKO.A = A_sKO;
    model_sKO.rhs = rhs_sKO;
    model_sKO.sense = sen_sKO;
    model_sKO.lb = lb_sKO;
    model_sKO.ub = ub_sKO;
    model_sKO.vtype = vtype_sKO;
    model_sKO.modelsense = 'min';  

    % --- Simulation 1: Maximize Biomass ---
    model_sKO_bio = model_sKO;
    model_sKO_bio.obj = [zeros(nRxns_sKO,1); kappa_val*ones(nRxns_sKO,1); kappa_val*ones(nRxns_sKO,1)];
    model_sKO_bio.obj(v0_sKO + biomass_idx) = -mu;  

    params.OutputFlag = 0;
    result_sKO_bio = gurobi(model_sKO_bio, params);

    if strcmp(result_sKO_bio.status, 'OPTIMAL')
        v_sol_sKO = result_sKO_bio.x(1:nRxns_sKO);
        alpha_sKO = result_sKO_bio.x(nRxns_sKO+1:2*nRxns_sKO);
        beta_sKO  = result_sKO_bio.x(2*nRxns_sKO+1:3*nRxns_sKO);

        growth_sKO_vbio(k) = v_sol_sKO(biomass_idx);
        prod_sKO_vbio(k) = v_sol_sKO(prod_idx);
        slack_sKO_vbio(k) = sum(alpha_sKO > 1e-6 | beta_sKO > 1e-6);

        if nKOs < nKOsLimit
            fprintf('  [Bio Max] → Biomass = %.4f, Product = %.4f, slacks = %d\n', ...
                        growth_sKO_vbio(k), prod_sKO_vbio(k), slack_sKO_vbio(k));
        end
    else
        if nKOs < nKOsLimit, fprintf('  [Bio Max] → Infeasible\n'); end
    end

    % --- Simulation 2: Maximize Product ---
    model_sKO_prod = model_sKO;
    model_sKO_prod.obj = [zeros(nRxns_sKO,1); kappa_val*ones(nRxns_sKO,1); kappa_val*ones(nRxns_sKO,1)];
    model_sKO_prod.obj(v0_sKO + prod_idx) = -mu;   

    result_sKO_prod = gurobi(model_sKO_prod, params);

    if strcmp(result_sKO_prod.status, 'OPTIMAL')
        v_sol_sKO = result_sKO_prod.x(1:nRxns_sKO);
        alpha_sKO = result_sKO_prod.x(nRxns_sKO+1:2*nRxns_sKO);
        beta_sKO  = result_sKO_prod.x(2*nRxns_sKO+1:3*nRxns_sKO);

        growth_sKO_prod(k) = v_sol_sKO(biomass_idx);
        prod_sKO_prod(k) = v_sol_sKO(prod_idx);
        slack_sKO_prod(k) = sum(alpha_sKO > 1e-6 | beta_sKO > 1e-6);

        if nKOs < nKOsLimit
            fprintf('  [Prd Max] → Biomass = %.4f, Product = %.4f, slacks = %d\n', ...
                        growth_sKO_prod(k), prod_sKO_prod(k), slack_sKO_prod(k));
        end
    else
        if nKOs < nKOsLimit, fprintf('  [Prd Max] → Infeasible\n'); end
    end
end

% Aggregate Final Results
sKO_results = table('Size', [length(KO_names), 5], ...
    'VariableTypes', {'cell', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'KO_names', 'Max_growth', 'Min_product', 'Min_growth', 'Max_product'});

sKO_results.KO_names    = KO_names;
sKO_results.Max_growth  = growth_sKO_vbio;
sKO_results.Min_product = prod_sKO_vbio;
sKO_results.Min_growth  = growth_sKO_prod;
sKO_results.Max_product = prod_sKO_prod;

exportFilename_sKO = strcat('Ecoli_strategies_sKO_', metName, '.xlsx');
% writetable(sKO_results, exportFilename_sKO, "FileType", "spreadsheet")

%% Analysis for Single-KO Strategies
fprintf('\n--- Starting single TF KO Analysis ---\n');
% Preallocate results table
sKO_results = table();
% Base WT reference for z
z0_base = z_struct.(best_threshold);

% Base PROM dimensions
[nMets_sKO, nRxns_sKO] = size(modelWT_FBA.S);
v0_sKO = 0; a0_sKO = nRxns_sKO; b0_sKO = 2 * nRxns_sKO;
nVars_sKO = 3 * nRxns_sKO;
vtype_sKO = repmat('C', nVars_sKO, 1);

for s = 1:height(resultsTableFilter)
    strat_name = sprintf('Multi_Iter%d', resultsTableFilter.Iteration(s));
    z_multi    = resultsTableFilter.z_sol{s};
    kappa_val  = resultsTableFilter.Kappa(s);
    orig_prod  = final_results.Max_product(s);
    
    % Find all perturbations in this strategy
    ko_indices  = find(z0_base == 1 & z_multi < 0.5);
    act_indices = find(z0_base == 0 & z_multi > 0.5);
    all_perturb_indices = [ko_indices; act_indices];
    
    if isempty(all_perturb_indices)
        continue;
    end
    fprintf('Analyzing %s (%d total perturbations)...\n', strat_name, length(all_perturb_indices));
    
    for p = 1:length(all_perturb_indices)
        tf_idx  = all_perturb_indices(p);
        tf_name = TF_list{tf_idx};
    
        % Single-TF simulation: Start from WT, apply ONLY this perturbation
        z_sKO = z0_base;   
        if z_multi(tf_idx) < 0.5
            z_sKO(tf_idx) = 0;
            perturb_type = 'Single_KO';
        else
            z_sKO(tf_idx) = 1;
            perturb_type = 'Single_Activation';
        end
        
        % Build problem
        TF_activity_sKO = TF_exp .* z_sKO;
        y_sKO = beta_matrix * TF_activity_sKO;
        y_sKO = gene_to_rxn_map * y_sKO;
        y_sKO = y_sKO ./ tau;
        y_sKO(isnan(y_sKO)) = 1; y_sKO(~isfinite(y_sKO)) = 1;
        
        A_sKO = sparse(0, nVars_sKO); rhs_sKO = []; sen_sKO = '';
        A_sKO = [A_sKO; sparse(modelWT_FBA.S), sparse(nMets_sKO, 2*nRxns_sKO)];
        rhs_sKO = [rhs_sKO; zeros(nMets_sKO,1)];
        sen_sKO = [sen_sKO; repmat('=', nMets_sKO, 1)];
        
        for i = 1:nRxns_sKO
            y_i = y_sKO(i);
            if any(gene_to_rxn_map(i,:)) && y_i < 0.99 
                row_lb = zeros(1, nVars_sKO); row_lb(v0_sKO + i) = 1; row_lb(a0_sKO + i) = 1;
                A_sKO = [A_sKO; row_lb]; rhs_sKO = [rhs_sKO; modelWT_FBA.lb(i) * y_i]; sen_sKO = [sen_sKO; '>'];
                row_ub = zeros(1, nVars_sKO); row_ub(v0_sKO + i) = 1; row_ub(b0_sKO + i) = -1;
                A_sKO = [A_sKO; row_ub]; rhs_sKO = [rhs_sKO; modelWT_FBA.ub(i) * y_i]; sen_sKO = [sen_sKO; '<'];
            end
        end  
        
        lb_sKO = [modelWT_FBA.lb; zeros(2*nRxns_sKO,1)]; 
        ub_sKO = [modelWT_FBA.ub; inf(2*nRxns_sKO,1)];   
        lb_sKO(v0_sKO + biomass_idx) = biomassWT_FBA * 0.5;
        ub_sKO(v0_sKO + biomass_idx) = 1000;
        
        % Set up optimization to Maximize Product
        model_sKO = struct('A', A_sKO, 'rhs', rhs_sKO, 'sense', sen_sKO, ...
                           'lb', lb_sKO, 'ub', ub_sKO, 'vtype', vtype_sKO, 'modelsense', 'max');
        
        % Reward product heavily to overcome kappa slack penalties
        weight_prod = max(1e4, kappa_val * 10);
        model_sKO.obj = [zeros(nRxns_sKO,1); -kappa_val*ones(nRxns_sKO,1); -kappa_val*ones(nRxns_sKO,1)];
        model_sKO.obj(v0_sKO + prod_idx) = weight_prod;   
        
        params.OutputFlag = 0;
        result_sKO = gurobi(model_sKO, params);
        
        if strcmp(result_sKO.status, 'OPTIMAL')
            v_sol_sKO = result_sKO.x(1:nRxns_sKO);
            sKO_prod = v_sol_sKO(prod_idx);
        else
            sKO_prod = NaN;
        end
        
        if ~isnan(sKO_prod) && orig_prod > 1e-6
            prod_drop_pct = 100 * (orig_prod - sKO_prod) / orig_prod;
        else
            prod_drop_pct = 0;
        end
        
        if prod_drop_pct > 5 
            gene_status = 'Essential Driver';
        elseif prod_drop_pct < -5 
            gene_status = 'Harmful (Collateral)';
        else
            gene_status = 'Redundant (Silent)';
        end
        
        sKO_results = [sKO_results; ...
            table({strat_name}, {tf_name}, {perturb_type}, orig_prod, sKO_prod, prod_drop_pct, {gene_status}, ...
            'VariableNames', {'Strategy', 'Perturbed_TF', 'Action', 'Multi_Product', 'sKO_Product', 'Product_Drop_Pct', 'Classification'})];
            
        fprintf('  [%s] %s: Prod went from %.4f -> %.4f (Drop: %5.1f%%) [%s]\n', ...
            perturb_type, tf_name, orig_prod, sKO_prod, prod_drop_pct, gene_status);
    end
end

fileName_sKO = strcat('Ecoli_singleKO_', metName{1,1}, '.xlsx');
% writetable(sKO_results, fileName_sKO, "FileType", "spreadsheet");

%% Leave-One-Out (LOO) Analysis for Multi-KO Strategies
fprintf('\n--- Starting Leave-One-Out (LOO) Analysis ---\n');
loo_results = table();
z0_base = z_struct.(best_threshold);

[nMets_loo, nRxns_loo] = size(modelWT_FBA.S);
v0_loo = 0; a0_loo = nRxns_loo; b0_loo = 2 * nRxns_loo;
nVars_loo = 3 * nRxns_loo;
vtype_loo = repmat('C', nVars_loo, 1);

for s = 1:height(resultsTableFilter)
    strat_name = sprintf('Multi_Iter%d', resultsTableFilter.Iteration(s));
    z_multi    = resultsTableFilter.z_sol{s};
    kappa_val  = resultsTableFilter.Kappa(s);
    orig_prod  = final_results.Max_product(s);
    
    ko_indices  = find(z0_base == 1 & z_multi < 0.5);
    act_indices = find(z0_base == 0 & z_multi > 0.5);
    all_perturb_indices = [ko_indices; act_indices];
    
    if isempty(all_perturb_indices), continue; end
    fprintf('Analyzing %s (%d total perturbations)...\n', strat_name, length(all_perturb_indices));
    
    for p = 1:length(all_perturb_indices)
        tf_idx = all_perturb_indices(p);
        tf_name = TF_list{tf_idx};
        
        % "Leave One Out" by restoring it to WT state
        z_loo = z_multi;
        z_loo(tf_idx) = z0_base(tf_idx); 
        
        perturb_type = 'KO_Restored';
        if z0_base(tf_idx) == 0, perturb_type = 'Activation_Restored'; end
        
        TF_activity_loo = TF_exp .* z_loo;
        y_loo = beta_matrix * TF_activity_loo;
        y_loo = gene_to_rxn_map * y_loo;
        y_loo = y_loo ./ tau;
        y_loo(isnan(y_loo)) = 1; y_loo(~isfinite(y_loo)) = 1;
       
        A_loo = sparse(0, nVars_loo); rhs_loo = []; sen_loo = '';
        A_loo = [A_loo; sparse(modelWT_FBA.S), sparse(nMets_loo, 2*nRxns_loo)];
        rhs_loo = [rhs_loo; zeros(nMets_loo,1)];
        sen_loo = [sen_loo; repmat('=', nMets_loo, 1)];
        
        for i = 1:nRxns_loo
            y_i = y_loo(i);
            if any(gene_to_rxn_map(i,:)) && y_i < 0.99 
                row_lb = zeros(1, nVars_loo); row_lb(v0_loo + i) = 1; row_lb(a0_loo + i) = 1;
                A_loo = [A_loo; row_lb]; rhs_loo = [rhs_loo; modelWT_FBA.lb(i) * y_i]; sen_loo = [sen_loo; '>'];
                row_ub = zeros(1, nVars_loo); row_ub(v0_loo + i) = 1; row_ub(b0_loo + i) = -1;
                A_loo = [A_loo; row_ub]; rhs_loo = [rhs_loo; modelWT_FBA.ub(i) * y_i]; sen_loo = [sen_loo; '<'];
            end
        end  
        
        lb_loo = [modelWT_FBA.lb; zeros(2*nRxns_loo,1)]; 
        ub_loo = [modelWT_FBA.ub; inf(2*nRxns_loo,1)];   
        lb_loo(v0_loo + biomass_idx) = biomassWT_FBA * 0.5;
        ub_loo(v0_loo + biomass_idx) = 1000;
        
        model_loo = struct('A', A_loo, 'rhs', rhs_loo, 'sense', sen_loo, ...
                           'lb', lb_loo, 'ub', ub_loo, 'vtype', vtype_loo, 'modelsense', 'max');
                       
        % Reward product heavily to overcome kappa slack penalties
        weight_prod = max(1e4, kappa_val * 10);
        model_loo.obj = [zeros(nRxns_loo,1); -kappa_val*ones(nRxns_loo,1); -kappa_val*ones(nRxns_loo,1)];
        model_loo.obj(v0_loo + prod_idx) = weight_prod;   
        
        params.OutputFlag = 0;
        result_loo = gurobi(model_loo, params);
        
        if strcmp(result_loo.status, 'OPTIMAL')
            v_sol_loo = result_loo.x(1:nRxns_loo);
            loo_prod = v_sol_loo(prod_idx);
        else
            loo_prod = NaN;
        end
        
        if ~isnan(loo_prod) && orig_prod > 1e-6
            prod_drop_pct = 100 * (orig_prod - loo_prod) / orig_prod;
        else
            prod_drop_pct = 0;
        end
        
        if prod_drop_pct > 5 
            gene_status = 'Essential Driver';
        elseif prod_drop_pct < -5 
            gene_status = 'Harmful (Collateral)';
        else
            gene_status = 'Redundant (Silent)';
        end
        
        loo_results = [loo_results; ...
            table({strat_name}, {tf_name}, {perturb_type}, orig_prod, loo_prod, prod_drop_pct, {gene_status}, ...
            'VariableNames', {'Strategy', 'Restored_TF', 'Action', 'Multi_Product', 'LOO_Product', 'Product_Drop_Pct', 'Classification'})];
            
        fprintf('  Restored %s: Prod went from %.4f -> %.4f (Drop: %5.1f%%) [%s]\n', ...
            tf_name, orig_prod, loo_prod, prod_drop_pct, gene_status);
    end
end

fileName_loo = strcat('Ecoli_leaveOneOut_', metName{1,1}, '.xlsx');
% writetable(loo_results, fileName_loo, "FileType", "spreadsheet");

%% Production envelopes analysis
fprintf('\n--- Starting Production Envelope Analysis ---\n');

% Resolution for simulations
n_steps_env = 20;  % 20 steps for production envelope
biomass_max = biomassWT_FBA; 

% Storage structure
data_multi  = struct('name', {}, 'env', {});
data_single = struct('name', {}, 'env', {});
data_loo    = struct('name', {}, 'env', {});

strategy_list = {};
base_kappa = WT_results.(best_threshold).kappa;
base_z     = z_struct.(best_threshold);

% 1. Add Multi KOs
for i = 1:height(resultsTableFilter)
    s_name = sprintf('Multi_Iter%d', resultsTableFilter.Iteration(i));
    s_z    = resultsTableFilter.z_sol{i};
    s_kappa = resultsTableFilter.Kappa(i);
    strategy_list(end+1,:) = {s_name, 'Multi', s_z, s_kappa};
end

% 2. Add Single KOs (ensure we only simulate unique single perturbations)
if height(sKO_results) > 0
    unique_sKO = unique(sKO_results(:, {'Perturbed_TF', 'Action'}));
    for i = 1:height(unique_sKO)
        tf_name = unique_sKO.Perturbed_TF{i};
        action  = unique_sKO.Action{i};
        tf_idx  = find(strcmp(TF_list, tf_name));
        
        if ~isempty(tf_idx)
            s_z = base_z;
            if contains(action, 'KO')
                s_z(tf_idx) = 0;
            else
                s_z(tf_idx) = 1;
            end
            s_name = sprintf('Single_%s', tf_name);
            strategy_list(end+1,:) = {s_name, 'Single', s_z, base_kappa};
        end
    end
end

% 3. Add Leave-One-Out (LOO) strategies
if height(loo_results) > 0
    for i = 1:height(loo_results)
        strat_name = loo_results.Strategy{i};
        tf_name    = loo_results.Restored_TF{i};
        
        % Find the original z_multi array
        iter_str = strrep(strat_name, 'Multi_Iter', '');
        iter_num = str2double(iter_str);
        multi_idx = find(resultsTableFilter.Iteration == iter_num, 1);
        
        if ~isempty(multi_idx)
            z_multi   = resultsTableFilter.z_sol{multi_idx};
            kappa_val = resultsTableFilter.Kappa(multi_idx);
            tf_idx    = find(strcmp(TF_list, tf_name));
            
            s_z = z_multi;
            s_z(tf_idx) = base_z(tf_idx); % Restore TF to WT condition
            
            s_name = sprintf('LOO_%s_%s', strat_name, tf_name);
            strategy_list(end+1,:) = {s_name, 'LOO', s_z, kappa_val};
        end
    end
end

% Simulation loop 
for s = 1:size(strategy_list,1)
    s_name   = strategy_list{s,1};
    s_type   = strategy_list{s,2};
    z_curr   = strategy_list{s,3};
    kap_curr = strategy_list{s,4};
    
    % Recalculate y based on the strategy's z configuration
    TF_activity_env = TF_exp .* z_curr;
    y_curr = beta_matrix * TF_activity_env;
    y_curr = gene_to_rxn_map * y_curr;
    y_curr = y_curr ./ tau;
    y_curr(isnan(y_curr)) = 1; 
    y_curr(~isfinite(y_curr)) = 1;
    
    [nMets_p, nRxns_p] = size(modelWT_FBA.S);
    nVars_p = 3 * nRxns_p; 
    v0_p = 0; a0_p = nRxns_p; b0_p = 2*nRxns_p;
    
    A_p = sparse(0, nVars_p); rhs_p = []; sen_p = '';
    A_p = [A_p; sparse(modelWT_FBA.S), sparse(nMets_p, 2*nRxns_p)];
    rhs_p = [rhs_p; zeros(nMets_p, 1)];
    sen_p = [sen_p; repmat('=', nMets_p, 1)];
    
    % Scale bounds using y
    for r = 1:nRxns_p
        y_val = y_curr(r);
        if any(gene_to_rxn_map(r,:)) && y_val < 0.99
            % LB Constraint
            row = zeros(1, nVars_p); row(v0_p + r) = 1; row(a0_p + r) = 1;
            A_p = [A_p; row]; rhs_p = [rhs_p; modelWT_FBA.lb(r) * y_val]; sen_p = [sen_p; '>'];
            % UB Constraint
            row = zeros(1, nVars_p); row(v0_p + r) = 1; row(b0_p + r) = -1;
            A_p = [A_p; row]; rhs_p = [rhs_p; modelWT_FBA.ub(r) * y_val]; sen_p = [sen_p; '<'];
        end
    end
    
    lb_p = [modelWT_FBA.lb; zeros(2*nRxns_p, 1)];
    ub_p = [modelWT_FBA.ub; inf(2*nRxns_p, 1)];
    vtype_p = repmat('C', nVars_p, 1);
    
    % Base model constraints for the envelope
    m_base = struct('A', A_p, 'rhs', rhs_p, 'sense', sen_p, 'lb', lb_p, 'ub', ub_p, ...
                    'vtype', vtype_p, 'modelsense', 'max');
                    
    % Weight the product heavily to overcome negative kappa slack penalties
    weight_prod = max(1e4, kap_curr * 10);
    m_base.obj = [zeros(nRxns_p,1); -kap_curr*ones(nRxns_p,1); -kap_curr*ones(nRxns_p,1)];
    m_base.obj(v0_p + prod_idx) = weight_prod;
    
    % Production envelope run
    env_points = zeros(n_steps_env, 2);
    for b_step = 1:n_steps_env
        bio_target = (b_step-1) * (biomass_max / (n_steps_env-1));
        m_env = m_base;
        
        % Force biomass to specific step on the envelope
        m_env.lb(v0_p + biomass_idx) = bio_target;
        m_env.ub(v0_p + biomass_idx) = bio_target;
        
        params.OutputFlag = 0;
        res_env = gurobi(m_env, params);
        
        if strcmp(res_env.status, 'OPTIMAL')
            env_points(b_step, :) = [bio_target, res_env.x(prod_idx)];
        else
            env_points(b_step, :) = [bio_target, NaN]; % Infeasible at this growth rate
        end
    end
    
    % Store results based on category
    if strcmp(s_type, 'Multi')
        idx = length(data_multi) + 1;
        data_multi(idx).name = s_name;
        data_multi(idx).env  = env_points;
    elseif strcmp(s_type, 'Single')
        idx = length(data_single) + 1;
        data_single(idx).name = s_name;
        data_single(idx).env  = env_points;
    elseif strcmp(s_type, 'LOO')
        idx = length(data_loo) + 1;
        data_loo(idx).name = s_name;
        data_loo(idx).env  = env_points;
    end
end

% --- Refactor results to tables inline ---

% 1. Multi KOs
env_rows_multi = [];
for s = 1:length(data_multi)
    strat_name = string(data_multi(s).name);
    env_data   = data_multi(s).env;
    if isempty(env_data)
        continue;
    end
    n_points = size(env_data,1);
    temp_table = table( ...
        repmat(strat_name, n_points, 1), ...
        env_data(:,1), ...
        env_data(:,2), ...
        'VariableNames', {'Strategy', 'Biomass', 'Product'});
    env_rows_multi = [env_rows_multi; temp_table];
end
table_env_multi = env_rows_multi;

% 2. Single KOs
env_rows_single = [];
for s = 1:length(data_single)
    strat_name = string(data_single(s).name);
    env_data   = data_single(s).env;
    if isempty(env_data)
        continue;
    end
    n_points = size(env_data,1);
    temp_table = table( ...
        repmat(strat_name, n_points, 1), ...
        env_data(:,1), ...
        env_data(:,2), ...
        'VariableNames', {'Gene', 'Biomass', 'Product'});
    env_rows_single = [env_rows_single; temp_table];
end
table_env_single = env_rows_single;

% 3. Leave-One-Out
env_rows_loo = [];
for s = 1:length(data_loo)
    strat_name = string(data_loo(s).name);
    env_data   = data_loo(s).env;
    if isempty(env_data)
        continue;
    end
    n_points = size(env_data,1);
    temp_table = table( ...
        repmat(strat_name, n_points, 1), ...
        env_data(:,1), ...
        env_data(:,2), ...
        'VariableNames', {'LOO_Perturbation', 'Biomass', 'Product'});
    env_rows_loo = [env_rows_loo; temp_table];
end
table_env_loo = env_rows_loo;

% Finalize and Export
metName = model.metNames(find(model.S(:,prod_idx)));

if ~isempty(table_env_multi)
    table_env_multi.Metabolite(:,1) = metName;
    fileName_prod_mKO = strcat('Prod_envelope_multi_', metName{1,1}, '.xlsx');
    % writetable(table_env_multi, fileName_prod_mKO, "FileType", "spreadsheet")
end

if ~isempty(table_env_single)
    table_env_single.Metabolite(:,1) = metName;
    fileName_prod_sKO = strcat('Prod_envelope_sKO_', metName{1,1}, '.xlsx');
    % writetable(table_env_single, fileName_prod_sKO, "FileType", "spreadsheet")
end

if ~isempty(table_env_loo)
    table_env_loo.Metabolite(:,1) = metName;
    fileName_prod_loo = strcat('Prod_envelope_LOO_', metName{1,1}, '.xlsx');
    % writetable(table_env_loo, fileName_prod_loo, "FileType", "spreadsheet")
end

%%
toc