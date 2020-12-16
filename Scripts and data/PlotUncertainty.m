function [] = PlotUncertainty(doFig56, doFig4, res)
close all
%PLOTUNCERTAINTY Summary of this function goes here
%   Detailed explanation goes here
if nargin<1, doFig56=0; end
if nargin<2, doFig4=0; end
if nargin<3, res=0.01; end
if ~doFig56
    baseFolder = './Results, Fig2-3' ;
    load('./Results, Fig2-3/simple,  opt-eSS(181.58).mat')
elseif  doFig56
    baseFolder = './Results, Fig5-6' ;
    load('./Results, Fig5-6/simple, opt-eSS(213.71).mat')
end

files=dir([baseFolder '/minmax/*.mat']);
modelName='LipolysisModel';
[model,data, lb, ub, nParams, expInd, stimulus, dgf] = Init(modelName, doFig56, 1);

stimulusHighRes=table();
ins=log10(unique([data.InVitro.FFA.Ins; data.InVitro.Glycerol.Ins]));
stimulusHighRes.Ins=[0 10.^(ins(2):res:ins(end)) 0]';
stimulusHighRes.Iso=[10*ones(height(stimulusHighRes)-1,1); 0];

if length(optParam)==length(lb)
    optParam=[optParam 1-eps];
end
bestParam=optParam;
cost = costfunction(bestParam,model, expInd,  data, stimulus, doFig56, 1);
fprintf('Total cost: %.2f, chi2: %.2f. Pass: %d\n',cost, chi2inv(0.95,dgf), cost<chi2inv(0.95,dgf))

allParams=[];
for i = fliplr(1:length(files))
    load([files(i).folder '/' files(i).name],'optParam');
    if length(optParam)==length(lb)
        optParam=[optParam 1-eps];
    end
    allParams(i,:)=optParam;
end
allParams=unique(allParams,'rows');

if doFig4
  disp('Running batch 1/4')
end
[DR, DRHSL, DRDiabetes, TS] = SimulateAllParams(bestParam, allParams, model, data, expInd, stimulusHighRes, doFig56);

if doFig4
  disp('Running batch 2/4')
  [DR_exclude1] = SimulateAllParams(bestParam, allParams, model, data, expInd, stimulusHighRes,doFig56, 'exclude1');
  disp('Running batch 3/4')
  
  [DR_exclude2] = SimulateAllParams(bestParam, allParams, model, data, expInd, stimulusHighRes,doFig56, 'exclude2');
  disp('Running batch 4/4')
  [~, ~, ~, TS_exclude3] = SimulateAllParams(bestParam, allParams, model, data, expInd, stimulusHighRes(1,:),doFig56, 'exclude3');
end


%% Do the plotting 
%Setup simulations and data (in vitro) 
close all
allInVitroData=struct();
allInVitroData.Normal=data.InVitro;
allInVitroData.Diabetes=data.InVitro_diabetes;

% Plot estimation simulation  and data
if ~doFig4
  PlotInVivo(data, TS.Gly.Normal, 2,'Glycerol')
  PlotInVitro(DR, allInVitroData, {'Normal'}, 2)
end
if contains(baseFolder, 'Fig2-3') % Plot validation simulation and data

  if ~doFig4
    PlotInVitro(DRHSL, allInVitroData, {'Normal'}, 3)
    exportgraphics(figure(2),sprintf('./%s/Fig2-Estimation.pdf',baseFolder))
    exportgraphics(figure(3),sprintf('./%s/Fig3-Validation.pdf',baseFolder))
  elseif doFig4
    DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso', 'Glycerol'}))=[];
    PlotInVitro(DR, allInVitroData, {'Normal'}, 4,[], 1)
    PlotInVitro(DR_exclude1, allInVitroData, {'Normal'}, 4,[], 3)
    PlotInVitro(DR_exclude2, allInVitroData, {'Normal'}, 4,[], 5)
    
    PlotInVivo(data, TS.Gly.Normal, 4,'Glycerol', {'Fig3Epi','Fig3EpiPhe'}, 2)
    PlotInVivo(data, TS_exclude3.Gly.Normal, 4,'Glycerol', {'Fig3Epi','Fig3EpiPhe'}, 4)
  end
    
    
elseif contains(baseFolder, 'Fig5-6') % Plot diabetes in vitro predictions.
    exportgraphics(figure(2),sprintf('./%s/FigS1-Estimation.pdf',baseFolder))
 
    DRDiabetes.Normal(:,~ismember(DRDiabetes.Normal.Properties.VariableNames,{'Ins','Iso','FFA','HSL'}))=[];
    DRDiabetes.Diabetes(:,~ismember(DRDiabetes.Diabetes.Properties.VariableNames,{'Ins','Iso','FFA', 'HSL'}))=[];
    DRDiabetes.Diabetes.HSL=nan(size(DRDiabetes.Diabetes.HSL));
    PlotInVitro(DRDiabetes, allInVitroData, {'Normal','Diabetes'}, 5)
    subplot(2,2,2)
    axsLine=findobj(gca,'type','line');
    legend([axsLine(8) axsLine(5)], {'non-diabetic','diabetic'}, 'location','best')
    exportgraphics(figure(5),sprintf('./%s/Fig5-Estimation-new.pdf',baseFolder))
    
    for i = 2:size(TS.FFA.Normal,2)
        maxPred=max(max([TS.FFA.Normal{:,i}; TS.FFA.Diabetes{:,i}]));
        TS.FFA.Normal{:,i}=TS.FFA.Normal{:,i}./maxPred*100;
        TS.FFA.Diabetes{:,i}=TS.FFA.Diabetes{:,i}./maxPred*100;
    end
    TS.FFA.Diabetes.Properties.VariableNames=strcat(TS.FFA.Diabetes.Properties.VariableNames,'Diab');
    TS.FFA=[TS.FFA.Normal, TS.FFA.Diabetes(:,2:end)];
    PlotInVivo(data.InVivo.Fig1.Time, TS.FFA, 6, 'FFA',{'Fig1', 'Fig2Epi', 'Fig1Diab', 'Fig2EpiDiab'})
    subplot(2,2,1)
    h=findobj(gca,'Type','line');
    h(7).YData=nan(size(h(7).YData));
    subplot(2,2,2)
    h=findobj(gca,'Type','line');
    h(7).YData=nan(size(h(7).YData));
    exportgraphics(figure(6),sprintf('./%s/Fig6-InVivo-FFA.pdf',baseFolder))
    
end
end

function [DR, DRHSL, DRDiabetes, TS] = SimulateAllParams(optParam, allParams, model, data, expInd, stimulusHighRes, doDiabetes, exclude)
if nargin < 8, exclude=''; end
switch exclude
  case 'exclude1', allParams(:,28)=inf; optParam(28)=inf;
  case 'exclude2', allParams(:,29)=inf; optParam(29)=inf;
  case 'exclude3', allParams(:,27)=inf; optParam(27)=inf;
end


params=[exp(optParam(1:expInd)) optParam(expInd+1:end-1)];
diab=optParam(end);

[DRbest]=simulateDR(model, params, expInd, doDiabetes, diab, stimulusHighRes);
DRtmp.high=DRbest.Normal;
DRtmp.low=DRbest.Normal;
if ~doDiabetes
  DRbest.Diabetes=DRbest.Normal; 
end

DRbest.Diabetes{:,3:end}=nan;
DRtmpD.high=DRbest.Diabetes;
DRtmpD.low=DRbest.Diabetes;

[TSBest] = SimulateInVivo(params, model, data.InVivo);
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

for i = 1:size(allParams,1)
    optParam=allParams(i,:);
    
    params=[exp(optParam(1:expInd)) optParam(expInd+1:end-1)];
    diab=optParam(end);
    Tmp=simulateDR(model, params, expInd, doDiabetes,  diab, stimulusHighRes);
    DRtmp.low{:,:}=min(DRtmp.low{:,:}, Tmp.Normal{:,:});
    DRtmp.high{:,:}=max(DRtmp.high{:,:}, Tmp.Normal{:,:});
    if doDiabetes
        DRtmpD.low{:,:}=min(DRtmpD.low{:,:}, Tmp.Diabetes{:,:});
        DRtmpD.high{:,:}=max(DRtmpD.high{:,:}, Tmp.Diabetes{:,:});
        
        TmpD = SimulateInVivo(params, model, data.InVivo, diab);
        TStmpD.Gly.low{:,:}=min(TStmpD.Gly.low{:,:}, TmpD.Gly{:,:});
        TStmpD.Gly.high{:,:}=max(TStmpD.Gly.high{:,:}, TmpD.Gly{:,:});
        TStmpD.FFA.low{:,:}=min(TStmpD.FFA.low{:,:}, TmpD.FFA{:,:});
        TStmpD.FFA.high{:,:}=max(TStmpD.FFA.high{:,:}, TmpD.FFA{:,:});
    end
    
    Tmp = SimulateInVivo(params, model, data.InVivo);
    TStmp.Gly.low{:,:}=min(TStmp.Gly.low{:,:}, Tmp.Gly{:,:});
    TStmp.Gly.high{:,:}=max(TStmp.Gly.high{:,:}, Tmp.Gly{:,:});
    TStmp.FFA.low{:,:}=min(TStmp.FFA.low{:,:}, Tmp.FFA{:,:});
    TStmp.FFA.high{:,:}=max(TStmp.FFA.high{:,:}, Tmp.FFA{:,:});
    
    fprintf('%i of %i \n',i,size(allParams,1))
end


DR.Normal=ConcatenateTableColums(DRbest.Normal, DRtmp.low(:,3:end));
DR.Normal=ConcatenateTableColums(DR.Normal, DRtmp.high(:,3:end));
DR.Diabetes=ConcatenateTableColums(DRbest.Diabetes, DRtmpD.low(:,3:end));
DR.Diabetes=ConcatenateTableColums(DR.Diabetes, DRtmpD.high(:,3:end));
DRHSL=DR;
DRHSL.Normal(:,~ismember(DRHSL.Normal.Properties.VariableNames,{'Ins','Iso','HSL'}))=[];
DRDiabetes=DR;

%Setup simulations and data (in vivo)
TS.Gly.Normal=ConcatenateTableColums(TSBest.Gly, TStmp.Gly.low(:,2:end));
TS.Gly.Normal=ConcatenateTableColums(TS.Gly.Normal, TStmp.Gly.high(:,2:end));

if isempty(exclude)
  TS.FFA.Normal=ConcatenateTableColums(TSBest.FFA, TStmp.FFA.low(:,2:end));
  TS.FFA.Normal=ConcatenateTableColums(TS.FFA.Normal, TStmp.FFA.high(:,2:end));
  
  TS.FFA.Diabetes=ConcatenateTableColums(TSBestD.FFA, TStmpD.FFA.low(:,2:end));
  TS.FFA.Diabetes=ConcatenateTableColums(TS.FFA.Diabetes, TStmpD.FFA.high(:,2:end));
  DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso','FFA', 'Glycerol','PKB'}))=[];
else
  DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso', 'Glycerol'}))=[];
end
end