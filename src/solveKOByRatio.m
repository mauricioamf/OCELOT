function [biomass, flux] = solveKOByRatio(modelKO, biomass_idx, WT_kappa, mu, gene_to_rxn_map, y_KO_table, TF, protectedRxnIdx)
    S_KO = sparse(modelKO.S);
    [m_KO, n_KO] = size(S_KO);
    lb_base = modelKO.lb;
    ub_base = modelKO.ub;

    y_rxn_KO = y_KO_table.(TF);
    y_WT = y_KO_table{:,'WT'};

    lbp = lb_base;
    ubp = ub_base;

    isProtected = false(n_KO,1);
    isProtected(protectedRxnIdx) = true;

    for i = 1:n_KO
        if isProtected(i)
            continue;
        end
        if ~any(gene_to_rxn_map(i,:))
            continue;
        end

        y_wt = y_WT(i);
        y_ko = y_rxn_KO(i);

        if ~isfinite(y_wt) || ~isfinite(y_ko) || abs(y_wt) <= 1e-9
            continue;
        end

        y_ratio = y_ko / y_wt;
        y_ratio = max(0, min(y_ratio, 1));

        if abs(lb_base(i)) > 1e-6 && abs(ub_base(i)) > 1e-6
            lbp(i) = max(lb_base(i) * y_ratio, -1000);
            ubp(i) = min(ub_base(i) * y_ratio,  1000);
        elseif lb_base(i) < -1e-6
            lbp(i) = max(lb_base(i) * y_ratio, -1000);
        elseif ub_base(i) > 1e-6
            ubp(i) = min(ub_base(i) * y_ratio,  1000);
        end
    end

    lbp(~isfinite(lbp)) = lb_base(~isfinite(lbp));
    ubp(~isfinite(ubp)) = ub_base(~isfinite(ubp));

    infeasible_idx = find(lbp > ubp);
    for k = infeasible_idx'
        lbp(k) = min(lbp(k), ubp(k) - 1e-6);
        ubp(k) = max(lbp(k), ubp(k) + 1e-6);
    end

    gp.A = sparse([ ...
        S_KO,         sparse(m_KO,n_KO),   sparse(m_KO,n_KO); ...
       -speye(n_KO), -speye(n_KO),         sparse(n_KO,n_KO); ...
        speye(n_KO),  sparse(n_KO,n_KO),  -speye(n_KO)]);

    gp.rhs        = [zeros(m_KO,1); -lbp; ubp];
    gp.sense      = [repmat('=',m_KO,1); repmat('<',n_KO,1); repmat('<',n_KO,1)];
    gp.modelsense = 'min';
    gp.obj        = [zeros(n_KO,1); WT_kappa*ones(n_KO,1); WT_kappa*ones(n_KO,1)];
    gp.obj(biomass_idx) = -mu;
    gp.lb         = [lb_base; zeros(2*n_KO,1)];
    gp.ub         = [ub_base; inf(2*n_KO,1)];

    gparams = struct();
    gparams.OutputFlag = 0;
    result_KO = gurobi(gp, gparams);

    if isfield(result_KO, 'x')
        flux = result_KO.x(1:n_KO);
        biomass = flux(biomass_idx);
    else
        flux = [];
        biomass = 0;
    end
end