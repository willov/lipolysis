function [cost] = costfunction(param_in,model, expInd, data, stimulus, doDiabetes, doPlot)

if (nargin < 6), doDiabetes = 0; end
if (nargin < 7), doPlot = 0; end

if iscolumn(param_in); param_in=param_in'; end

params=[exp(param_in(1:expInd)) param_in(expInd+1:end)];
if length(param_in)==33 % hardcorded for LipolysisModel model, to account for the fact that some parameter sets where collected without diabetes as a free parameter. TODO: Extend previously collected parameter sets
    diabetes=params(end);
    params(end)=[];
else
    diabetes=0.66; % not used
end

costInVitro = CostInVivo(params, model, data, doPlot);
costTot = CostDR(model, params, expInd, stimulus, data,diabetes, doDiabetes, doPlot);

cost=costInVitro+costTot;

end



