function upsilon = buildUpsilon(gene_to_rxn_map, beta_prime, y_WT_base, protectedRxnIdx, upsilonMode)
    upsilon = gene_to_rxn_map * beta_prime;

    switch lower(string(upsilonMode))
        case "global_max"
            upsilon = upsilon ./ max(upsilon,[],'all');
            % disp("GLOBAL_MAX")
        otherwise
            for i = 1:size(upsilon,1)
                if y_WT_base(i) > 1e-9
                    upsilon(i,:) = upsilon(i,:) ./ y_WT_base(i);
                else
                    upsilon(i,:) = 0;
                end
            end
            % disp("WT_RATIO")
    end

    % upsilon = max(upsilon, 0);
    % upsilon = min(upsilon, 1);
    % 
    if ~isempty(protectedRxnIdx)
        upsilon(protectedRxnIdx,:) = 0;
    end
end