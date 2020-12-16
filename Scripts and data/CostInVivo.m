function [cost, responses] = CostInVivo(params, model, data, printCost)

try
ind = ismember(data.InVivo.Fig1.Time(1):data.InVivo.Fig1.Time(30), data.InVivo.Fig1.Time);

[responses] = SimulateInVivo(params, model, data.InVivo);
yh=responses.Gly.Fig1(ind);
y=data.InVivo.Fig1.Mean;
SEM=data.InVivo.Fig1.SEM;
costFig1 = sum((yh-y).^2 ./SEM.^2);

yh=[responses.Gly.Fig2Epi(ind); responses.Gly.Fig2Iso(ind)];
y=[data.InVivo.Fig2Epi.Mean; data.InVivo.Fig2Iso.Mean];
SEM=[data.InVivo.Fig2Epi.SEM; data.InVivo.Fig2Iso.SEM];
costFig2 = sum((yh-y).^2 ./SEM.^2);

yh=[responses.Gly.Fig3Epi(ind); responses.Gly.Fig3EpiPhe(ind)];
y=[data.InVivo.Fig3Epi.Mean; data.InVivo.Fig3EpiPhe.Mean];
SEM=[data.InVivo.Fig3Epi.SEM; data.InVivo.Fig3EpiPhe.SEM];
costFig3 = sum((yh-y).^2 ./SEM.^2);

cost=costFig1+costFig2+costFig3;
    
catch err
    cost=1e91; %returns a "large" error.
    responses=[];
end

if printCost
    fprintf('Fig1: %.4f\n', costFig1)
    fprintf('Fig2: %.4f\n', costFig2)
    fprintf('Fig3: %.4f\n', costFig3)
end
end

