function vec = replaceInfWithFinite(vec)
    if isempty(vec)
        return;
    end
    finiteVals = vec(isfinite(vec));
    if isempty(finiteVals)
        vec(:) = 0;
        return;
    end
    replacement = min(finiteVals);
    vec(isinf(vec)) = replacement;
    vec(isnan(vec)) = 0;
end