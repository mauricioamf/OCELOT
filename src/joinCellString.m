function out = joinCellString(c)
    if isempty(c)
        out = "";
    else
        out = string(strjoin(cellstr(string(c)), ', '));
    end
end