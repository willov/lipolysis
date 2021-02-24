function [] = PlotUncertainty(doDiabetes,res, modelName)
clear mex

if nargin<1, doDiabetes=0; end
if nargin<2, res=0.01; end
if nargin<3, modelName='lipolysis'; end

if ~doDiabetes
    baseFolder = './Parameter sets' ;
elseif  doDiabetes
    baseFolder = './Parameter sets (with diabetes)' ;
end

files=dir(sprintf('%s/PPL/*%s*.mat', baseFolder, modelName));

[model,data, lb, ub, nParams, expInd, stimulus, dgf] = Init(modelName, doDiabetes, 0);
limit  = chi2inv(0.95,dgf);

stimulusHighRes=table();
ins=log10(unique([data.InVitro.FFA.Ins; data.InVitro.Glycerol.Ins]));
stimulusHighRes.Ins=[0 10.^(ins(2):res:ins(end)) 0]';
stimulusHighRes.Iso=[0.01*ones(height(stimulusHighRes)-1,1); 0]; %10 nM = 0.01 µM

load(FindBestParametersFile(baseFolder, 1, [modelName ', opt-eSS']), 'optParam')

if length(optParam)==length(IQMparameters(model))-5
    optParam=[optParam 0];
end

cost = costfunction(optParam,model, expInd,  data, stimulus, doDiabetes);
fprintf('Total cost: %.2f, chi2: %.2f. Pass: %d\n',cost, limit, cost<limit)
params=[exp(optParam(1:expInd)) optParam(expInd+1:end)];

diab=params(end);
params(end)=[];

Best=simulateInVitro(model, params, expInd, diab, stimulusHighRes, 0);
InVitrotmp.high=Best.Normal;
InVitrotmp.low=Best.Normal;

if doDiabetes
    InVitrotmpD.high=Best.Diabetes;
    InVitrotmpD.low=Best.Diabetes;
end

TSBest = SimulateInVivo(params, model, data.InVivo, 0, 0);
TSBest.FFA{:,2:end}=TSBest.FFA{:,2:end}./TSBest.FFA{1,2:end}; %ta  bort?
TStmp.Gly.low=TSBest.Gly;
TStmp.Gly.high=TSBest.Gly;
TStmp.FFA.low=TSBest.FFA;
TStmp.FFA.high=TSBest.FFA;

TSBestD=TSBest;
TSBestD.Gly{:,2:end}=nan;
TSBestD.FFA{:,2:end}=nan;
TStmpD.Gly.low=TSBestD.Gly;
TStmpD.Gly.high=TSBestD.Gly;
TStmpD.FFA.low=TSBestD.FFA;
TStmpD.FFA.high=TSBestD.FFA;

optParams=[];
for i = fliplr(1:length(files))
    load([files(i).folder '/' files(i).name],'optParam');
    if length(optParam)==length(lb)
        optParam=[optParam 0];
    end
    optParams(i,:)=optParam;
end
optParams=unique(optParams,'rows');
if size(optParams,1)>0
    disp('Simulating parameter sets')
end
for i = 1:size(optParams,1)
    optParam=optParams(i,:);
    cost = costfunction(optParam,model, expInd,  data, stimulus, doDiabetes);
    if cost<limit+0.1
        params=[exp(optParam(1:expInd)) optParam(expInd+1:end-1)];
        diab=optParam(end);
        
        [Tmp]=simulateInVitro(model, params, expInd, diab, stimulusHighRes);
        InVitrotmp.low{:,:}=min(InVitrotmp.low{:,:}, Tmp.Normal{:,:});
        InVitrotmp.high{:,:}=max(InVitrotmp.high{:,:}, Tmp.Normal{:,:});
        
        TmpTS = SimulateInVivo(params, model, data.InVivo, 0, 0);
        if isfield(Tmp,'Diabetes')
            InVitrotmpD.low{:,:}=min(InVitrotmpD.low{:,:}, Tmp.Diabetes{:,:});
            InVitrotmpD.high{:,:}=max(InVitrotmpD.high{:,:}, Tmp.Diabetes{:,:});
            
            TmpTSD = SimulateInVivo(params, model, data.InVivo, diab, 0);
            TmpTSD.FFA{:,2:end}=TmpTSD.FFA{:,2:end}./TmpTS.FFA{1,2:end};
            
            TStmpD.Gly.low{:,:}=min(TStmpD.Gly.low{:,:}, TmpTSD.Gly{:,:});
            TStmpD.Gly.high{:,:}=max(TStmpD.Gly.high{:,:}, TmpTSD.Gly{:,:});
            TStmpD.FFA.low{:,:}=min(TStmpD.FFA.low{:,:}, TmpTSD.FFA{:,:});
            TStmpD.FFA.high{:,:}=max(TStmpD.FFA.high{:,:}, TmpTSD.FFA{:,:});
        end
        TmpTS.FFA{:,2:end}=TmpTS.FFA{:,2:end}./TmpTS.FFA{1,2:end};
        TStmp.Gly.low{:,:}=min(TStmp.Gly.low{:,:}, TmpTS.Gly{:,:});
        TStmp.Gly.high{:,:}=max(TStmp.Gly.high{:,:}, TmpTS.Gly{:,:});
        TStmp.FFA.low{:,:}=min(TStmp.FFA.low{:,:}, TmpTS.FFA{:,:});
        TStmp.FFA.high{:,:}=max(TStmp.FFA.high{:,:}, TmpTS.FFA{:,:});
        
    else
        disp("Not a valid solution")
    end
    if i==1
        fprintf('%i of %i \n|',i,size(optParams,1))
    elseif mod(i,25)==0
        fprintf(' %i of %i \n',i,size(optParams,1))
    else
        fprintf('|')
    end
end

%% Do the plotting
%Setup simulations and data (in vitro)
allInVitroData=struct();
allInVitroData.Normal=data.InVitro;

InVitro.Normal=ConcatenateTableColums(Best.Normal, InVitrotmp.low(:,3:end));
InVitro.Normal=ConcatenateTableColums(InVitro.Normal, InVitrotmp.high(:,3:end));
InVitroHSL=InVitro;

if doDiabetes
    allInVitroData.Diabetes=data.InVitro_diabetes;
    
    InVitro.Diabetes=ConcatenateTableColums(Best.Diabetes, InVitrotmpD.low(:,3:end));
    InVitro.Diabetes=ConcatenateTableColums(InVitro.Diabetes, InVitrotmpD.high(:,3:end));
    InVitroDiabetes=InVitro;
    
    TS.FFA.Diabetes=ConcatenateTableColums(TSBestD.FFA, TStmpD.FFA.low(:,2:end));
    TS.FFA.Diabetes=ConcatenateTableColums(TS.FFA.Diabetes, TStmpD.FFA.high(:,2:end));
end

%Setup simulations and data (in vivo)
TS.Gly.Normal=ConcatenateTableColums(TSBest.Gly, TStmp.Gly.low(:,2:end));
TS.Gly.Normal=ConcatenateTableColums(TS.Gly.Normal, TStmp.Gly.high(:,2:end));

TS.FFA.Normal=ConcatenateTableColums(TSBest.FFA, TStmp.FFA.low(:,2:end));
TS.FFA.Normal=ConcatenateTableColums(TS.FFA.Normal, TStmp.FFA.high(:,2:end));

InVitro.Normal(:,~ismember(InVitro.Normal.Properties.VariableNames,{'Ins','Iso','FFA', 'Glycerol','PKB'}))=[];
Best.Normal(:,~ismember(Best.Normal.Properties.VariableNames,{'Ins','Iso','FFA', 'Glycerol','PKB'}))=[];

%% Plot
if contains(modelName,'_noIns')
    PlotInVivo(data, TS.Gly.Normal, 52,'Glycerol')
    PlotInVitro(InVitro, allInVitroData, {'Normal'}, 52)
elseif ~contains(baseFolder, 'with diabetes') % Plot validation simulation and data
    
    PlotInVivo(data, TS.Gly.Normal, 51,'Glycerol')
    PlotInVitro(InVitro, allInVitroData, {'Normal'}, 51)
    
    data.InVivo.Fig3Epi=[];
    TS.Gly.Normal.Fig3Epi=[];
    PlotInVivo(data, TS.Gly.Normal, 3,'Glycerol')
    PlotInVitro(InVitro, allInVitroData, {'Normal'}, 3)
    
    InVitroHSL.Normal(:,~ismember(InVitroHSL.Normal.Properties.VariableNames,{'Ins','Iso','HSL'}))=[];
    PlotInVitro(InVitroHSL, allInVitroData, {'Normal'}, 4)
elseif contains(baseFolder, 'with diabetes') % Plot diabetes in vitro predictions.
    
    PlotInVivo(data, TS.Gly.Normal, 53,'Glycerol')
    PlotInVitro(InVitro, allInVitroData, {'Normal'}, 53)
    
    data.InVivo.Fig3Epi=[];
    TS.Gly.Normal.Fig3Epi=[];
    PlotInVivo(data, TS.Gly.Normal, 54,'Glycerol')
    PlotInVitro(InVitro, allInVitroData, {'Normal'}, 54)
    
    
    InVitroDiabetes.Normal(:,~ismember(InVitroDiabetes.Normal.Properties.VariableNames,{'Ins','Iso','FFA', 'Glycerol', 'HSL', 'Reesterification'}))=[];
    InVitroDiabetes.Diabetes(:,~ismember(InVitroDiabetes.Diabetes.Properties.VariableNames,{'Ins','Iso','FFA', 'Glycerol', 'HSL', 'Reesterification'}))=[];
    
    allInVitroData.Normal.cAMP=[];
    PlotInVitro(InVitroDiabetes, allInVitroData, {'Normal','Diabetes'}, 6)
    
    %%
    if height(stimulusHighRes)<800
        fprintf('\n\nNote: running with a low resolution might yield slightly incorrect bounds for the reesterification\n')
    end
    reestMax = max(max(InVitroDiabetes.Normal.Reesterification(1:end-1,:)));
    reestMin = min(min(InVitroDiabetes.Normal.Reesterification(1:end-1,:)));
    fprintf('Normal reesterification (iso+ins stimulated) interval: %.2f - %.2f\n',reestMin, reestMax)
    
    reestDiabMax = max(max(InVitroDiabetes.Diabetes.Reesterification(1:end-1,:)));
    reestDiabMin = min(min(InVitroDiabetes.Diabetes.Reesterification(1:end-1,:)));
    fprintf('Diabetic reesterification (iso+ins stimulated) interval: %.2f - %.2f\n\n', reestDiabMin, reestDiabMax)
    
    %%
    TS.FFA.Diabetes.Properties.VariableNames=strcat(TS.FFA.Diabetes.Properties.VariableNames,'Diab');
    TS.FFA=[TS.FFA.Normal, TS.FFA.Diabetes(:,2:end)];
    PlotInVivo(data.InVivo.Fig1.Time, TS.FFA, 7, 'FFA',{'Fig1', 'Fig2Epi', 'Fig1Diab', 'Fig2EpiDiab'})
    subplot(2,2,1)
    ylabel({'Released FA', 'fold over non-diabetic basal'})
    subplot(2,2,2)
    ylabel({'Released FA', 'fold over non-diabetic basal'})
    
end
end

