function [modelAuto, biomass_idx, biomassWT_FBA, FBAsol, protectedRxnIdx, protectedRxns] = configureAutotrophicModel(model, params)
    modelAuto = model;

    biomass_idx = find(strcmp(modelAuto.rxns, params.biomassRxn), 1);
    if isempty(biomass_idx)
        error('Biomass reaction "%s" not found in the model.', string(params.biomassRxn));
    end

    modelAuto = changeObjective(modelAuto, char(params.biomassRxn));

    offBiomass = toStringArray(params.autotrophicOffBiomassRxns);
    for i = 1:numel(offBiomass)
        if any(strcmp(modelAuto.rxns, offBiomass(i)))
            modelAuto = changeRxnBounds(modelAuto, char(offBiomass(i)), 0, 'b');
        end
    end

    if ~isempty(params.photonRxn) && any(strcmp(modelAuto.rxns, params.photonRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.photonRxn), -1000, 'l');
    end
    if ~isempty(params.glucoseRxn) && any(strcmp(modelAuto.rxns, params.glucoseRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.glucoseRxn), 0, 'b');
    end
    if ~isempty(params.hco3Rxn) && any(strcmp(modelAuto.rxns, params.hco3Rxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.hco3Rxn), params.HCO3Uptake, 'l');
    end
    if ~isempty(params.co2Rxn) && any(strcmp(modelAuto.rxns, params.co2Rxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.co2Rxn), 0, 'b');
    end
    if ~isempty(params.h2co3TransportRxn) && any(strcmp(modelAuto.rxns, params.h2co3TransportRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.h2co3TransportRxn), 0, 'u');
    end

    disableRxns = toStringArray(params.autotrophicDisableRxns);
    for i = 1:numel(disableRxns)
        if any(strcmp(modelAuto.rxns, disableRxns(i)))
            modelAuto = changeRxnBounds(modelAuto, char(disableRxns(i)), 0, 'b');
        end
    end

    if ~isempty(params.fnorRxn) && any(strcmp(modelAuto.rxns, params.fnorRxn))
        modelAuto = changeRxnBounds(modelAuto, char(params.fnorRxn), 0, 'l');
    end

    modelAuto = applyBoundChangeList(modelAuto, params.additionalBoundChanges);

    solAuto1 = optimizeCbModel(modelAuto, 'max');
    if isempty(solAuto1) || ~isfield(solAuto1, 'f') || isempty(solAuto1.f)
        error('Autotrophic WT optimization failed at the first biomass step.');
    end

    modelAuto2 = changeObjective(modelAuto, char(params.photonRxn), 1);
    modelAuto2 = changeRxnBounds(modelAuto2, char(params.biomassRxn), solAuto1.f, 'b');
    solAuto2 = optimizeCbModel(modelAuto2, 'max');

    if isempty(solAuto2) || ~isfield(solAuto2, 'f') || isempty(solAuto2.f)
        error('Autotrophic photon minimization failed.');
    end

    modelAuto = changeRxnBounds(modelAuto, char(params.photonRxn), solAuto2.f, 'l');
    solAuto3 = optimizeCbModel(modelAuto, 'max');

    if isempty(solAuto3) || ~isfield(solAuto3, 'x') || isempty(solAuto3.x)
        error('Autotrophic WT optimization failed after fixing photon uptake.');
    end

    biomassWT_FBA = solAuto3.x(biomass_idx);
    FBAsol = solAuto3;

    protectedRxns = unique([ ...
        toStringArray(params.protectedRxns); ...
        toStringArray(params.autotrophicOffBiomassRxns); ...
        string(params.photonRxn); string(params.glucoseRxn); string(params.hco3Rxn); ...
        string(params.co2Rxn); string(params.h2co3TransportRxn); string(params.fnorRxn); ...
        toStringArray(params.autotrophicDisableRxns)]);

    protectedRxns = protectedRxns(strlength(protectedRxns) > 0);
    protectedRxnIdx = find(ismember(string(modelAuto.rxns), protectedRxns));
end