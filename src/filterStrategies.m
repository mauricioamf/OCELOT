function resultsTableFilter = filterStrategies(resultsTable)
    if isempty(resultsTable)
        resultsTableFilter = table();
        return;
    end

    numericVars = vartype('numeric');
    numericData = table2array(resultsTable(:, numericVars));
    rowHasNaN = any(isnan(numericData), 2);
    rowNoPerturb = (resultsTable.NumKO == 0) & (resultsTable.NumActivated == 0);
    rowAllKOorAllAct = (resultsTable.NumKO == numel(resultsTable.z0{1})) | (resultsTable.NumActivated == numel(resultsTable.z0{1}));
    rowZallZeroes = cellfun(@(v) isempty(v) || all(v == 0), resultsTable.z_sol);

    badRows = rowHasNaN | rowNoPerturb | rowAllKOorAllAct | rowZallZeroes;
    fprintf('Filtering %d / %d rows (NaN / no-perturb / all-TF).\n', sum(badRows), height(resultsTable));
    resultsTableFilter = resultsTable(~badRows, :);
end