function [cost_tot, cost, responses]=CostDR(model, param, expInd, stimulus, data, diab, doDiabetes, printCost)
if nargin <7, doDiabetes=0; end 
try

    [responses,~, basal]=simulateDR(model, param, expInd, doDiabetes, diab, stimulus);
    
    if any(abs(basal{:,:})<1e-3)
        penalty=sum(1./abs(basal{:,:}));
    else
        penalty=0;
    end
    
    costGlycerol=0;
    costcAMP=0;
    costFFA=0;
    costHSL=0;
    costPKB=0;
    costFFADiabetes=0;
    
    residuals=(responses.Normal{:,'Glycerol'}-data.InVitro.Glycerol.Mean).^2./data.InVitro.Glycerol.SEM.^2;
    costGlycerol=sum(residuals(~any(isnan(data.InVitro.Glycerol{:,3:4}),2))); 
    
    residuals=(responses.Normal{:,'FFA'}-data.InVitro.FFA.Mean).^2./data.InVitro.FFA.SEM.^2;
    costFFA=sum(residuals(~any(isnan(data.InVitro.FFA{:,3:4}),2))); 
    
    ind=ismember(responses.Normal{:,{'Ins','Iso'}},data.InVitro.PKB{:,{'Ins','Iso'}},'rows');
    residuals=(responses.Normal{ind,'PKB'}-data.InVitro.PKB.Mean).^2./data.InVitro.PKB.SEM.^2;
    costPKB=sum(residuals(~any(isnan(data.InVitro.PKB{:,3:4}),2))); 
    
    if doDiabetes
        residuals=(responses.Normal{:,'HSL'}-data.InVitro.HSL.Mean).^2./data.InVitro.HSL.SEM.^2;
        costHSL=sum(residuals(~any(isnan(data.InVitro.HSL{:,3:4}),2))); 
        
        residuals=(responses.Diabetes{:,'FFA'}-data.InVitro_diabetes.FFA.Mean).^2./data.InVitro_diabetes.FFA.SEM.^2;
        costFFADiabetes=sum(residuals(~any(isnan(data.InVitro_diabetes.FFA{:,3:4}),2)));
    end
    
    cost=costGlycerol+costFFA+costPKB+costcAMP+costHSL+costFFADiabetes;
    cost_tot=cost+penalty;
    
    %% Costs below here are not used when estimating.
    residuals=(responses.Normal{:,'HSL'}-data.InVitro.HSL.Mean).^2./data.InVitro.HSL.SEM.^2;
    costHSL=sum(residuals(~any(isnan(data.InVitro.HSL{:,3:4}),2))); 
    
    residuals=(responses.Normal{:,'cAMP'}-data.InVitro.cAMP.Mean).^2./data.InVitro.cAMP.SEM.^2;
    costcAMP=sum(residuals(~any(isnan(data.InVitro.cAMP{:,3:4}),2))); 
    
    if  ~doDiabetes
      costFFADiabetes = nan;
    end
catch 
  cost=1e90; %returns a "large" error.
  cost_tot=1e91;
  costGlycerol=inf;
  costcAMP=inf;
  costFFA=inf;
    costHSL=inf;
    costPKB=inf;
    responses=[];
end

if printCost
    fprintf('DR, Glycerol: %.2f\n',costGlycerol)
    fprintf('DR, FFA: %.2f\n',costFFA)
    fprintf('DR, HSL: %.2f\n',costHSL)
    fprintf('DR, cAMP: %.2f\n',costcAMP)
    fprintf('DR, PKB: %.2f\n',costPKB)
    fprintf('DR, FFA-diabetes: %.2f\n',costFFADiabetes)

end

end

