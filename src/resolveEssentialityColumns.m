function [geneCol, labelCol] = resolveEssentialityColumns(essentialityTable, params)
    geneCol = '';
    labelCol = '';

    if isfield(params, 'essentialityGeneColumn') && strlength(string(params.essentialityGeneColumn)) > 0
        geneCol = char(params.essentialityGeneColumn);
    end
    if isfield(params, 'essentialityLabelColumn') && strlength(string(params.essentialityLabelColumn)) > 0
        labelCol = char(params.essentialityLabelColumn);
    end

    vars = essentialityTable.Properties.VariableNames;
    if isempty(geneCol)
        geneCandidates = {'Gene_ID','GeneID','Gene','gene','TF','tf','ORF','orf','locus_tag'};
        for i = 1:numel(geneCandidates)
            if any(strcmp(vars, geneCandidates{i}))
                geneCol = geneCandidates{i};
                break;
            end
        end
        if isempty(geneCol)
            for i = 1:numel(vars)
                col = essentialityTable.(vars{i});
                if iscellstr(col) || isstring(col) || ischar(col)
                    geneCol = vars{i};
                    break;
                end
            end
        end
    end

    if isempty(labelCol)
        labelCandidates = {'Essentiality_0_11notkeeping','Essentiality','essentiality','is_essential','label','Label'};
        for i = 1:numel(labelCandidates)
            if any(strcmp(vars, labelCandidates{i}))
                labelCol = labelCandidates{i};
                break;
            end
        end
        if isempty(labelCol)
            for i = 1:numel(vars)
                if strcmp(vars{i}, geneCol)
                    continue;
                end
                col = essentialityTable.(vars{i});
                if isnumeric(col) || islogical(col)
                    u = unique(col(~isnan(double(col))));
                    if numel(u) <= 3
                        labelCol = vars{i};
                        break;
                    end
                end
            end
        end
    end

    if isempty(geneCol) || isempty(labelCol)
        error('Could not resolve essentiality gene/label columns. Set params.essentialityGeneColumn and params.essentialityLabelColumn explicitly.');
    end
end