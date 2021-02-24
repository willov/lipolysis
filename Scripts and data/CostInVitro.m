function [cost, responses, costHSL]=CostInVitro(model, param, expInd, stimulus, data, diab, doDiabetes)
if nargin <7, doDiabetes=0; end
try
    
    [responses, steady, ~, basal]=simulateInVitro(model, param, expInd, diab, stimulus);
    
    if any(abs(basal{:,:})<1e-3)
        penalty=sum(1./abs(basal{:,:}))+50;
    else
        penalty = 0;
    end
    
    costHSL=0;
    costReesterificationDiab = 0;
    
    residuals=(responses.Normal{:,'Glycerol'}-data.InVitro.Glycerol.Mean).^2./data.InVitro.Glycerol.SEM.^2;
    costGlycerol=sum(residuals(~any(isnan(data.InVitro.Glycerol{:,3:4}),2))); %only som residuals where
    
    residuals=(responses.Normal{:,'FFA'}-data.InVitro.FFA.Mean).^2./data.InVitro.FFA.SEM.^2;
    costFFA=sum(residuals(~any(isnan(data.InVitro.FFA{:,3:4}),2))); %only som residuals where
    
    ind=ismember(responses.Normal{:,{'Ins'}},data.InVitro.PKB{:,{'Ins'}},'rows');
    residuals=(responses.Normal{ind,'PKB'}-data.InVitro.PKB.Mean).^2./data.InVitro.PKB.SEM.^2;
    costPKB=sum(residuals(~any(isnan(data.InVitro.PKB{:,3:4}),2))); %only som residuals where%
    
    if doDiabetes
        residuals=(responses.Normal{:,'HSL'}-data.InVitro.HSL.Mean).^2./data.InVitro.HSL.SEM.^2;
        costHSL=sum(residuals(~any(isnan(data.InVitro.HSL{:,3:4}),2))); %only som residuals where%
        
        ind=ismember(responses.Normal{:,{'Ins'}},data.InVitro_diabetes.Reesterification{:,{'Ins'}},'rows');
        residuals=(responses.Diabetes{ind,'Reesterification'}-data.InVitro_diabetes.Reesterification.Mean).^2./data.InVitro_diabetes.Reesterification.SEM.^2;
        costReesterificationDiab=sum(residuals(~any(isnan(data.InVitro_diabetes.Reesterification{:,3:4}),2))); %only som re
    end
    
    cost=costGlycerol+costFFA+costPKB+costHSL+costReesterificationDiab;
    cost=cost+penalty+steady;
    
catch err
    cost=1e90; %returns a "large" error.
    costHSL=inf;
    responses=[];
end



end

