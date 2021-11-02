function [cost, responses] = CostInVivo(params, model, data, tol, useFig2EpiData)
if nargin<4, tol = 0; end
try
    ind = ismember(data.Fig1.Time(1):data.Fig1.Time(30), data.Fig1.Time);
    
    [responses, steady] = SimulateInVivo(params, model, data, 0, tol);
    
    %% Fig1
    yh=responses.Gly.Fig1(ind);
    y=data.Fig1.Mean;
    SEM=data.Fig1.SEM;
    costFig1 = sum((yh-y).^2 ./SEM.^2);
    
    if useFig2EpiData
        %% Fig2
        yh=[responses.Gly.Fig2Epi(ind); responses.Gly.Fig2Iso(ind)];
        y=[data.Fig2Epi.Mean; data.Fig2Iso.Mean];
        SEM=[data.Fig2Epi.SEM; data.Fig2Iso.SEM];
        costFig2 = sum((yh-y).^2 ./SEM.^2);
        
        %% Fig3
        yh=responses.Gly.Fig3EpiPhe(ind);
        y=data.Fig3EpiPhe.Mean;
        SEM=data.Fig3EpiPhe.SEM;
        costFig3 = sum((yh-y).^2 ./SEM.^2);
    else
        %% Fig2
        yh=responses.Gly.Fig2Iso(ind);
        y=data.Fig2Iso.Mean;
        SEM=data.Fig2Iso.SEM;
        costFig2 = sum((yh-y).^2 ./SEM.^2);

        %% Fig3
        yh=[responses.Gly.Fig3Epi(ind); responses.Gly.Fig3EpiPhe(ind)];
        y=[data.Fig3Epi.Mean; data.Fig3EpiPhe.Mean];
        SEM=[data.Fig3Epi.SEM; data.Fig3EpiPhe.SEM];
        costFig3 = sum((yh-y).^2 ./SEM.^2);
    end

    %% Total cost
    cost=costFig1+costFig2+costFig3+steady;
    
catch err
    cost=1e91; %returns a "large" error.
    responses=[];
end


end

