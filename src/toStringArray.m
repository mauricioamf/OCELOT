function out = toStringArray(x)
    if isempty(x)
        out = strings(0,1);
    else
        out = string(x(:));
    end
    out = out(strlength(out) > 0);
end