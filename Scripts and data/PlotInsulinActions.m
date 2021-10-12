function [] = PlotInsulinActions(useFig2EpiData)
if nargin<1, useFig2EpiData=1; end

modelName = 'lipolysis';
if useFig2EpiData
    baseFolder='./Parameter sets';
else
    baseFolder='./Parameter sets (with alternative epi data)';
end

[model,data, lb, ub, nParams, expInd, stimulus, dgf] = Init(modelName, 0, 0);
inVitroData.Normal=data.InVitro;


stimulusHighRes=table();
ins=log10(unique([data.InVitro.FFA.Ins; data.InVitro.Glycerol.Ins]));
stimulusHighRes.Ins=[0 10.^(ins(2):.01:ins(end)) 0]';
stimulusHighRes.Iso=[0.01*ones(height(stimulusHighRes)-1,1); 0]; %10 nM = 0.01 µM

figure(5)
clf
%% Plot original model
load(FindBestParametersFile(baseFolder, 1, [modelName ', opt-eSS']), 'optParam')

diab=0;

params=[exp(optParam(1:expInd)) optParam(expInd+1:end)];
DR=simulateInVitro(model, params, expInd, diab, stimulusHighRes, 0);
TS = SimulateInVivo(params, model, data.InVivo, 0, 0);
DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso','Glycerol'}))=[];

PlotInVivo(data, TS.Gly, 5,'Glycerol', {'Fig3Epi','Fig3EpiPhe'},0, [2,2])
PlotInVitro(DR, inVitroData, {'Normal'}, 5, 1)

%% Plot without action 1)
params_no1 = params;
params_no1(strcmp(IQMparameters(model),'EC501'))=inf;

DR=simulateInVitro(model, params_no1, expInd, diab, stimulusHighRes, 0);
DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso','Glycerol'}))=[];
PlotInVitro(DR, inVitroData, {'Normal'}, 5, 3)

%% Plot without action 2)
params_no2 = params;
params_no2(strcmp(IQMparameters(model),'EC502'))=inf;

DR=simulateInVitro(model, params_no2, expInd, diab, stimulusHighRes, 0);
DR.Normal(:,~ismember(DR.Normal.Properties.VariableNames,{'Ins','Iso','Glycerol'}))=[];
PlotInVitro(DR, inVitroData, {'Normal'}, 5, 5)
%% Plot without action 3)
params_no3 = params;
params_no3(strcmp(IQMparameters(model),'EC503'))=inf;

TS = SimulateInVivo(params_no3, model, data.InVivo, 0, 0);
PlotInVivo(data, TS.Gly, 5,'Glycerol', {'Fig3Epi','Fig3EpiPhe'},0, [4,4])
end

