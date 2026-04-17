%% Running OCELOT for TF essentiality prediction

% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);

clear

load('../manuscript/Models/iML1515.mat');

expression_data = readtable('../manuscript/Data/Ecoli_count_matrix.csv', ...
    'ReadRowNames', true, 'VariableNamingRule', 'preserve');
beta_df = readtable('../manuscript/Data/Ecoli_GRN_lasso.csv', ...
    'ReadRowNames', true);
essentialityTable = readtable('../manuscript/Data/Ecoli_gene_essentiality.csv');

params.selectedReference = "control__wt_glc";
params.growthType = 'heterotrophic';
params.biomassRxn = 'BIOMASS_Ec_iML1515_core_75p37M';
params.glucoseRxn = 'EX_glc__D_e';
params.glucoseUptake = -10;
params.mediumRxns = ["EX_pi_e","EX_co2_e","EX_fe3_e","EX_h_e", ...
                     "EX_mn2_e","EX_fe2_e","EX_zn2_e","EX_mg2_e", ...
                     "EX_ca2_e","EX_ni2_e","EX_cu2_e","EX_sel_e", ...
                     "EX_cobalt2_e","EX_h2o_e","EX_mobd_e","EX_so4_e", ...
                     "EX_nh4_e","EX_k_e","EX_na1_e","EX_cl_e", ...
                     "EX_o2_e","EX_tungs_e","EX_slnt_e"];

params.performEssentialityEvaluation = true;
params.essentialityTable = essentialityTable;
params.essentialityGeneColumn = 'Gene_ID';
params.essentialityLabelColumn = 'Essentiality_0_1';
params.cutoffPerc = 0.1;
params.showPlots = "cutoff";

params.runKnockout = true;
params.runEngineering = false;

results = OCELOT(model, expression_data, beta_df, params);

%% Running OCELOT for predicting metabolic engineering strategies

% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP');
% changeCobraSolverParams('LP', 'feasTol', 1e-9);

clear

load('../manuscript/Models/iML1515.mat');

expression_data = readtable('../manuscript/Data/Ecoli_count_matrix.csv', ...
    'ReadRowNames', true, 'VariableNamingRule', 'preserve');
beta_df = readtable('../manuscript/Data/Ecoli_GRN_lasso.csv', ...
    'ReadRowNames', true);
essentialityTable = readtable('../manuscript/Data/Ecoli_gene_essentiality.csv');

params.selectedReference = "control__wt_glc";
params.growthType = 'heterotrophic';
params.biomassRxn = 'BIOMASS_Ec_iML1515_core_75p37M';
params.glucoseRxn = 'EX_glc__D_e';
params.glucoseUptake = -10;
params.mediumRxns = ["EX_pi_e","EX_co2_e","EX_fe3_e","EX_h_e", ...
                     "EX_mn2_e","EX_fe2_e","EX_zn2_e","EX_mg2_e", ...
                     "EX_ca2_e","EX_ni2_e","EX_cu2_e","EX_sel_e", ...
                     "EX_cobalt2_e","EX_h2o_e","EX_mobd_e","EX_so4_e", ...
                     "EX_nh4_e","EX_k_e","EX_na1_e","EX_cl_e", ...
                     "EX_o2_e","EX_tungs_e","EX_slnt_e"];

params.runKnockout = false;
params.runEngineering = true;
params.prodRxn = 'EX_ac_e';
params.maxPerturbations = 5;
params.maxIterations = 12;
params.flexBiomass = 0.5;
params.flexProduct = 0.5;
params.validateStrategies = true;
params.validateSingleKOs = true;
params.validateLeaveOneOut = true;

results_Eng = OCELOT(model, expression_data, beta_df, params);