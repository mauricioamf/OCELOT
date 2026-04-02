function [biomass, gp, result, kappa] = solveProblemWithScaledBounds(model, biomass_idx, yvec, gene_to_rxn_map, protectedRxnIdx, mu, kappaStart, kappaMin, kappaDecay)
    S = sparse(model.S);
    [m, n] = size(S);
    lb_base = model.lb;
    ub_base = model.ub;

    lbp = lb_base;
    ubp = ub_base;

    isProtected = false(n,1);
    isProtected(protectedRxnIdx) = true;

    yvec(~isfinite(yvec)) = 1;

    for i = 1:n
        if isProtected(i)
            continue;
        end
        if ~any(gene_to_rxn_map(i,:))
            continue;
        end

        y = yvec(i);

        if abs(lb_base(i)) > 1e-6 && abs(ub_base(i)) > 1e-6
            lbp(i) = max(lb_base(i) * y, -1000);
            ubp(i) = min(ub_base(i) * y,  1000);
        elseif lb_base(i) < -1e-6
            lbp(i) = max(lb_base(i) * y, -1000);
        elseif ub_base(i) > 1e-6
            ubp(i) = min(ub_base(i) * y,  1000);
        end
    end

    lbp(~isfinite(lbp)) = lb_base(~isfinite(lbp));
    ubp(~isfinite(ubp)) = ub_base(~isfinite(ubp));

    infeasible_idx = find(lbp > ubp);
    for k = infeasible_idx'
        lbp(k) = min(lbp(k), ubp(k) - 1e-6);
        ubp(k) = max(lbp(k), ubp(k) + 1e-6);
    end

    found_solution = false;
    gp = [];
    result = [];
    biomass = NaN;
    kappa = NaN;

    kappaLoop = kappaStart;
    while kappaLoop >= kappaMin
        gp.A = sparse([ ...
            S,          sparse(m,n),      sparse(m,n); ...
           -speye(n),  -speye(n),         sparse(n,n); ...
            speye(n),   sparse(n,n),     -speye(n)]);

        gp.rhs        = [zeros(m,1); -lbp; ubp];
        gp.sense      = [repmat('=',m,1); repmat('<',n,1); repmat('<',n,1)];
        gp.modelsense = 'min';
        gp.obj        = [zeros(n,1); kappaLoop*ones(n,1); kappaLoop*ones(n,1)];
        gp.obj(biomass_idx) = -mu;
        gp.lb         = [lb_base; zeros(2*n,1)];
        gp.ub         = [ub_base; inf(2*n,1)];

        gparams = struct();
        gparams.OutputFlag = 0;
        resultLoop = gurobi(gp, gparams);

        if isfield(resultLoop, 'x')
            v_sol = resultLoop.x(1:n);
            biomassLoop = v_sol(biomass_idx);

            if biomassLoop > 1e-6
                fprintf('→ predicted biomass = %.4f, kappa = %.3e\n', biomassLoop, kappaLoop);
                biomass = biomassLoop;
                result = resultLoop;
                kappa = kappaLoop;
                found_solution = true;
                break;
            end
        end

        kappaLoop = kappaLoop / kappaDecay;
    end

    if ~found_solution
        error('No feasible solution was found.');
    end
end