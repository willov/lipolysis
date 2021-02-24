function [cost] = costfunction(param_in,model, expInd, data, stimulus, doDiabetes)

if (nargin < 6), doDiabetes = 0; end

if iscolumn(param_in); param_in=param_in'; end

params=[exp(param_in(1:expInd)) param_in(expInd+1:end)];
if length(param_in)+4==length(IQMparameters(model))
    diabetes=params(end);
    params(end)=[];
else
    diabetes=0.66; % not used
end

costInVivo = CostInVivo(params, model, data.InVivo, 0);%% InVivo data
costInVitro = CostInVitro(model, params, expInd, stimulus, data,diabetes, doDiabetes);%% InVitro data

cost=costInVivo+costInVitro;

end
