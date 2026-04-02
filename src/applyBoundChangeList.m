function model = applyBoundChangeList(model, changes)
    if isempty(changes)
        return;
    end

    if ~isstruct(changes)
        error('additionalBoundChanges must be a struct array with fields rxn, value, and bound.');
    end

    for i = 1:numel(changes)
        if ~isfield(changes(i), 'rxn') || ~isfield(changes(i), 'value') || ~isfield(changes(i), 'bound')
            error('Each bound-change entry must contain fields rxn, value, and bound.');
        end
        if any(strcmp(model.rxns, changes(i).rxn))
            model = changeRxnBounds(model, char(changes(i).rxn), changes(i).value, changes(i).bound);
        end
    end
end