function [avg_expr_data, unique_conditions, selectedWTReference] = preprocessExpressionData(expression_data, selectedReference)
    all_samples = expression_data.Properties.VariableNames;
    condition_names = regexprep(all_samples, '__\d+$', '');

    for i = 1:numel(condition_names)
        if ~isempty(condition_names{i}) && isstrprop(condition_names{i}(1), 'digit')
            condition_names{i} = ['X' condition_names{i}];
        end
    end

    unique_conditions = unique(condition_names, 'stable');
    wt_idx = find(contains(unique_conditions, selectedReference), 1, 'first');
    if isempty(wt_idx)
        error('WT reference "%s" not found in expression_data.', string(selectedReference));
    end
    selectedWTReference = unique_conditions{wt_idx};

    avg_expr_data = array2table( ...
        zeros(size(expression_data,1), numel(unique_conditions)), ...
        'RowNames', expression_data.Properties.RowNames, ...
        'VariableNames', unique_conditions);

    for c = 1:numel(unique_conditions)
        cond = unique_conditions{c};
        idx = strcmp(condition_names, cond);
        avg_expr_data{:,c} = mean(expression_data{:,idx}, 2, 'omitnan');
    end
end