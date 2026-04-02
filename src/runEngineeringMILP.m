function [resultsTable, results_all] = runEngineeringMILP(modelWT_FBA, biomass_idx, prod_idx, biomassWT_FBA, maxProduct, TF_list, z0, upsilon, kappa, params)
    [nMets, nRxns] = size(modelWT_FBA.S);
    nTF = numel(TF_list);

    v0     = 0;
    a0     = v0 + nRxns;
    b0     = a0 + nRxns;
    z0_off = b0 + nRxns;
    nVars  = z0_off + nTF;

    vtype = [repmat('C', nRxns,1); repmat('C', nRxns,1); repmat('C', nRxns,1); repmat('B', nTF,1)];

    gparams = struct();
    gparams.OutputFlag = 0;
    gparams.TimeLimit = params.timeLimit;
    gparams.MIPGap = params.mipGap;
    gparams.MIPFocus = params.mipFocus;

    flexBiomass = params.flexBiomass(:)';
    flexProduct = params.flexProduct(:)';

    if numel(flexBiomass) == 1 && numel(flexProduct) > 1
        flexBiomass = repmat(flexBiomass, size(flexProduct));
    elseif numel(flexProduct) == 1 && numel(flexBiomass) > 1
        flexProduct = repmat(flexProduct, size(flexBiomass));
    elseif numel(flexBiomass) ~= numel(flexProduct)
        error('flexBiomass and flexProduct must have the same length, or one must be scalar.');
    end

    results_all = {};
    integer_cuts = {};

    for idxFlex = 1:numel(flexBiomass)
        fb = flexBiomass(idxFlex);
        fp = flexProduct(idxFlex);

        min_biomass_step3 = biomassWT_FBA * fb;
        min_product_step3 = maxProduct * fp;

        lb_base = modelWT_FBA.lb;
        ub_base = modelWT_FBA.ub;
        lb_base(biomass_idx) = min_biomass_step3;
        lb_base(prod_idx) = min_product_step3;

        A = sparse(0, nVars);
        rhs = [];
        sen = '';

        A = [A; sparse(modelWT_FBA.S), sparse(nMets, nVars - nRxns)];
        rhs = [rhs; zeros(nMets,1)];
        sen = [sen; repmat('=', nMets, 1)];

        for i = 1:nRxns
            if any(upsilon(i,:))
                row = zeros(1, nVars);
                row(v0 + i) = 1;
                row(a0 + i) = 1;
                row(z0_off+1 : z0_off+nTF) = -lb_base(i) * upsilon(i,:);
                A = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '>'];

                row = zeros(1, nVars);
                row(v0 + i) = 1;
                row(b0 + i) = -1;
                row(z0_off+1 : z0_off+nTF) = ub_base(i) * upsilon(i,:);
                A = [A; row];
                rhs = [rhs; 0];
                sen = [sen; '<'];
            end
        end

        perturb_row = zeros(1, nVars);
        for j = 1:nTF
            if z0(j) == 1
                perturb_row(z0_off + j) = -1;
            else
                perturb_row(z0_off + j) = 1;
            end
        end
        rhs_perturb = params.maxPerturbations - sum(z0 == 1);

        A = [A; perturb_row];
        rhs = [rhs; rhs_perturb];
        sen = [sen; '<'];

        lb = [lb_base; zeros(2*nRxns,1); zeros(nTF,1)];
        ub = [ub_base; inf(2*nRxns,1); ones(nTF,1)];

        for iteration = 1:params.maxIterations
            A_with_cuts = A;
            rhs_with_cuts = rhs;
            sen_with_cuts = sen;

            for cut_idx = 1:numel(integer_cuts)
                cut = integer_cuts{cut_idx};
                row = zeros(1, nVars);
                for j = 1:nTF
                    if cut.z_pattern(j) == 1
                        row(z0_off + j) = 1;
                    else
                        row(z0_off + j) = -1;
                    end
                end
                A_with_cuts = [A_with_cuts; row];
                rhs_with_cuts = [rhs_with_cuts; cut.sum_value - 2];
                sen_with_cuts = [sen_with_cuts; '<'];
            end

            obj = zeros(nVars,1);
            for j = 1:nTF
                if z0(j) == 1
                    obj(z0_off + j) = 1;
                else
                    obj(z0_off + j) = -1;
                end
            end
            obj(a0+1 : b0) = -kappa;
            obj(b0+1 : z0_off) = -kappa;

            model_step3 = struct();
            model_step3.A = sparse(A_with_cuts);
            model_step3.rhs = rhs_with_cuts;
            model_step3.sense = sen_with_cuts;
            model_step3.lb = lb;
            model_step3.ub = ub;
            model_step3.vtype = vtype;
            model_step3.obj = obj;
            model_step3.modelsense = 'max';

            result_step3 = gurobi(model_step3, gparams);

            if isfield(result_step3, 'status') && strcmp(result_step3.status, 'OPTIMAL')
                z_sol = result_step3.x(z0_off+1 : z0_off+nTF);
                v_sol = result_step3.x(1:nRxns);
                alpha_vals = result_step3.x(a0+1 : b0);
                beta_vals = result_step3.x(b0+1 : z0_off);

                knockout_idx = find(z0 == 1 & z_sol < 0.5);
                activated_idx = find(z0 == 0 & z_sol > 0.5);

                knocked_out_TFs = TF_list(knockout_idx);
                activated_TFs = TF_list(activated_idx);

                biomass_eng = v_sol(biomass_idx);
                product_eng = v_sol(prod_idx);
                num_slacks_used = sum(alpha_vals > 1e-6 | beta_vals > 1e-6);

                results_all(end+1,:) = { ...
                    fb, fp, biomass_eng, product_eng, ...
                    knocked_out_TFs, activated_TFs, ...
                    numel(knocked_out_TFs), numel(activated_TFs), ...
                    num_slacks_used, kappa, z0, z_sol, iteration, 'Original'}; %#ok<AGROW>

                cut_pattern = struct();
                cut_pattern.z_pattern = round(z_sol);
                cut_pattern.sum_value = sum(cut_pattern.z_pattern);

                is_new_pattern = true;
                for cut_idx = 1:numel(integer_cuts)
                    if isequal(integer_cuts{cut_idx}.z_pattern, cut_pattern.z_pattern)
                        is_new_pattern = false;
                        break;
                    end
                end
                if is_new_pattern
                    integer_cuts{end+1} = cut_pattern; %#ok<AGROW>
                end
            else
                if iteration == 1
                    results_all(end+1,:) = {fb, fp, NaN, NaN, {}, {}, 0, 0, NaN, kappa, z0, NaN, 0, 'Original'}; %#ok<AGROW>
                end
                break;
            end
        end
    end

    if isempty(results_all)
        resultsTable = table();
        return;
    end

    resultsTable = cell2table(results_all, ...
        'VariableNames', {'BiomassFlex','ProductFlex','Biomass','Product', ...
                          'Knockouts','Activated','NumKO','NumActivated', ...
                          'SlacksUsed','Kappa','z0','z_sol','Iteration','Strategy'});

    resultsTable = resultsTable(:, {'Strategy','Iteration','BiomassFlex','ProductFlex','Kappa','Biomass','Product', ...
                                    'Knockouts','Activated','NumKO','NumActivated','SlacksUsed','z0','z_sol'});
end