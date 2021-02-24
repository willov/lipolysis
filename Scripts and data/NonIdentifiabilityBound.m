function [bound] = NonIdentifiabilityBound(paramName,bound, direction)
if strcmp(direction,'min')
    if regexp(paramName, '^k')
        bound = log(exp(bound)/1e3);
    elseif contains(paramName, 'EC50')
        bound = -6;
    else
        bound = 0;
    end
else
    if regexp(paramName, '^k')
        bound = log(exp(bound)*1e3);
    elseif strcmp(paramName, 'phe_effect')
        bound = 1;
    elseif strcmp(paramName, 'isoscale')
        bound = 100;
    elseif contains(paramName, 'min')
        bound = 100;
    elseif contains(paramName, 'EC50')
        bound = 4;
    elseif ismember(paramName,{'n1','n2','n3'})
        bound = 4;
    end
end
end