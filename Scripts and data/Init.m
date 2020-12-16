function [model,data, lb, ub, nParams,expInd, stimulus, dgf] = Init(modelName, doDiabetes, doMEX)
if nargin<1; modelName='LipolysisModel'; end
if nargin<2, doDiabetes=0; end
if nargin<3, doMEX=0; end

if doMEX
    disp('Compiling the model, please wait')
    optModel = IQMmodel([modelName '.txt']);
    IQMmakeMEXmodel(optModel,modelName);
end
model=str2func(modelName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CREATE THE EXPDATA STRUCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./Data/stichfig1.mat'); stichfig1.ins = stichfig1.ins * 6.945*1e-3; % converting from mU/L to nM
load('./Data/stichfig2.mat'); stichfig2.ins = stichfig2.ins * 6.945*1e-3;
load('./Data/stichfig3.mat'); stichfig3.ins = stichfig3.ins * 6.945*1e-3;

data.InVivo.Fig1.Mean=stichfig1.val(1:30);
data.InVivo.Fig1.SEM=stichfig1.se(1:30);
data.InVivo.Fig1.Time=stichfig1.time(1:30);
data.InVivo.Fig1.Ins=stichfig1.ins;

data.InVivo.Fig2Epi.Mean=stichfig2.epi(1:30);
data.InVivo.Fig2Epi.SEM=stichfig2.epi_se(1:30);
data.InVivo.Fig2Epi.Time=stichfig2.time(1:30);
data.InVivo.Fig2Epi.Ins=stichfig2.ins(1:30);

data.InVivo.Fig2Iso.Mean=stichfig2.iso(1:30);
data.InVivo.Fig2Iso.SEM=stichfig2.iso_se(1:30);
data.InVivo.Fig2Iso.Time=stichfig2.time(1:30);
data.InVivo.Fig2Iso.Ins=stichfig2.ins(1:30);

data.InVivo.Fig3Epi.Mean=stichfig3.epi(1:30);
data.InVivo.Fig3Epi.SEM=stichfig3.epi_se(1:30);
data.InVivo.Fig3Epi.Time=stichfig3.time(1:30);
data.InVivo.Fig3Epi.Ins=stichfig3.ins(1:30);

data.InVivo.Fig3EpiPhe.Mean=stichfig3.epiphe(1:30);
data.InVivo.Fig3EpiPhe.SEM=stichfig3.epiphe_se(1:30);
data.InVivo.Fig3EpiPhe.Time=stichfig3.time(1:30);
data.InVivo.Fig3EpiPhe.Ins=stichfig3.ins(1:30);


load 'Data/expData'
expData.FFA.SEM(expData.FFA.SEM==0)=nan;
expData.Glycerol.SEM(expData.Glycerol.SEM==0)=nan;
expData.HSL.SEM(expData.HSL.SEM==0)=nan;
expData.cAMP.SEM(expData.cAMP.SEM==0)=nan;
expData.PKB.SEM(expData.PKB.SEM==0)=nan;

expData.FFA.Mean(2)=nan;
expData.FFA.SEM(2)=nan;

data.InVitro.FFA=expData.FFA;
data.InVitro.Glycerol=expData.Glycerol;
data.InVitro.cAMP=expData.cAMP;
data.InVitro.HSL=expData.HSL;
data.InVitro.PKB=expData.PKB;

%% Diabetes data
load 'Data/expDataDiabetes'
data.InVitro_diabetes.FFA=expDataDiabetes.FFA; 

%% DGF

dgfcAMP=sum(~isnan(data.InVitro.cAMP.SEM));
dgfGlycerol=sum(~isnan(data.InVitro.Glycerol.SEM));
dgfFFA=sum(~isnan(data.InVitro.FFA.SEM));
dgfHSL=sum(~isnan(data.InVitro.HSL.SEM));
dgfPKB=sum(~isnan(data.InVitro.PKB.SEM));
dgfFFADiabetes = sum(~isnan(data.InVitro_diabetes.FFA.SEM));

dgfFig1Epi=sum(~isnan(data.InVivo.Fig1.SEM));
dgfFig2Epi=sum(~isnan(data.InVivo.Fig2Epi.SEM));
dgfFig2Iso=sum(~isnan(data.InVivo.Fig2Iso.SEM));
dgfFig3Epi=sum(~isnan(data.InVivo.Fig3Epi.SEM));
dgfFig3Phe=sum(~isnan(data.InVivo.Fig3EpiPhe.SEM));


dgfInVitro=dgfGlycerol+dgfFFA+dgfPKB;%+dgfcAMP
if doDiabetes
    dgfInVitro=dgfInVitro+dgfFFADiabetes+dgfHSL;
end
dgfInVivo=dgfFig1Epi+dgfFig2Epi+dgfFig2Iso+dgfFig3Epi+dgfFig3Phe;

dgf=dgfInVivo+dgfInVitro-1-3; %-1 for drift parameter, -3 for scaling, -1 for experimental scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          THE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pNames, ~] = IQMparameters(model);
expInd=sum(~cellfun(@isempty,(regexp(pNames,'^k.+'))));
lb=repmat(log(1e-5), 1, expInd); % set lower bound for all "InVitro" parameters
ub=-1*lb; % set upper bound for all "InVitro" paramters


lb=[lb 0.1   0   0   0       -10   log10(0.5) -10 0  0 0]; %p=0?
ub=[ub 100 100 100 100        10   log10(1.1)  10 10 10 10];
nParams=length(lb);

%% setup stimulus
stimulus=table();
ins=unique([expData.FFA.Ins; expData.Glycerol.Ins]);
stimulus.Ins=[ins; 0];
stimulus.Iso=[10*ones(height(stimulus)-1,1); 0];

end

