function [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureBaseModel(model, params)
    growthType = lower(string(params.growthType));

    switch growthType
        case "autotrophic"
            [modelWT_FBA, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureAutotrophicModel(model, params);
        otherwise
            modelWT_FBA = model;

            exchangeIndex = find(contains(lower(string(modelWT_FBA.rxnNames)), "exchange"));
            exchangeIDs = string(modelWT_FBA.rxns(exchangeIndex));
            if params.removeLastExchange && ~isempty(exchangeIDs)
                exchangeIDs(end) = [];
            end

            biomass_idx = find(strcmp(modelWT_FBA.rxns, params.biomassRxn), 1);
            if isempty(biomass_idx)
                error('Biomass reaction "%s" not found in the model.', string(params.biomassRxn));
            end

            if ~isempty(exchangeIDs)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, cellstr(exchangeIDs), 0, 'l');
            end
            if ~isempty(params.mediumRxns)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, cellstr(toStringArray(params.mediumRxns)), -1000, 'l');
            end
            modelWT_FBA = changeObjective(modelWT_FBA, modelWT_FBA.rxns{biomass_idx});

            if ~isempty(params.glucoseRxn)
                modelWT_FBA = changeRxnBounds(modelWT_FBA, char(params.glucoseRxn), params.glucoseUptake, 'b');
            end

            modelWT_FBA = applyBoundChangeList(modelWT_FBA, params.additionalBoundChanges);

            FBAsol = optimizeCbModel(modelWT_FBA);
            biomassWT_FBA = FBAsol.f;

            protectedRxns = unique([toStringArray(params.protectedRxns); strings(0,1)]);
            protectedRxnIdx = find(ismember(string(modelWT_FBA.rxns), protectedRxns));
    end
end