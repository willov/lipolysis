function [model,data, lb, ub, nParams,expInd, stimulus, dgf, opts] = Init(modelName, doDiabetes, doMEX)
if nargin<1; modelName='lipolysis'; end
if nargin<2, doDiabetes=0; end
if nargin<3, doMEX=0; end
format short e

if doMEX
    disp('Compiling the model, please wait')
    optModel = IQMmodel([modelName '.txt']);
    IQMmakeMEXmodel(optModel,modelName);
end
model=str2func(modelName);

%%  Setup in vivo data
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

%% Setup in vitro data
load('Data/expData_fig14','expData')
expData.FFA.SEM(end-2:end-1)= mean(expData.FFA.SEM([2:end-3 end]));

expData.FFA.SEM(expData.FFA.SEM==0)=nan;
expData.Glycerol.SEM(expData.Glycerol.SEM==0)=nan;
expData.HSL.SEM(expData.HSL.SEM==0)=nan;
expData.PKB.SEM(expData.PKB.SEM==0)=nan;

data.InVitro.FFA=expData.FFA;
data.InVitro.Glycerol=expData.Glycerol;
data.InVitro.HSL=expData.HSL;
data.InVitro.PKB=expData.PKB;

%% Setup diabetes data
load('Data/expDataDiabetes','expDataDiabetes')
data.InVitro_diabetes.FFA=expDataDiabetes.FFA;
data.InVitro_diabetes.Glycerol=expDataDiabetes.Glycerol;
data.InVitro_diabetes.HSL=expDataDiabetes.HSL;
data.InVitro_diabetes.Reesterification=expDataDiabetes.reesterificationDiabetes;
data.InVitro.Reesterification=expDataDiabetes.reesterification;

%% setup degrees of freedom
dgfGlycerol=sum(~isnan(data.InVitro.Glycerol.SEM));
dgfFFA=sum(~isnan(data.InVitro.FFA.SEM));
dgfHSL=sum(~isnan(data.InVitro.HSL.SEM));
dgfPKB=sum(~isnan(data.InVitro.PKB.SEM));
dgfReesterificationDiab=sum(~isnan(data.InVitro_diabetes.Reesterification.SEM));

dgfFig1Epi=sum(~isnan(data.InVivo.Fig1.SEM));
dgfFig2Epi=sum(~isnan(data.InVivo.Fig2Epi.SEM));
dgfFig2Iso=sum(~isnan(data.InVivo.Fig2Iso.SEM));
dgfFig3Phe=sum(~isnan(data.InVivo.Fig3EpiPhe.SEM));

dgfInVitro=dgfGlycerol+dgfFFA+dgfPKB;
if doDiabetes
    dgfInVitro=dgfInVitro+dgfHSL+dgfReesterificationDiab;
end
dgfInVivo=dgfFig1Epi+dgfFig2Iso+dgfFig3Phe+dgfFig2Epi;

dgf=dgfInVivo+dgfInVitro-1-6; %-1 for drift parameter, -6 for scaling

%% Parameter bound setup
[pNames, ~] = IQMparameters(model);
expInd=sum(~cellfun(@isempty,(regexp(pNames,'^k.+'))));
lb=repmat(log(1e-6), 1, expInd); % set lower bound for all "InVitro" parameters
ub=repmat(log(1e6), 1, expInd); % set upper bound for all "InVitro" paramters

lb=[lb 0.6   8   0   0   0   -5   log(0.5) -5 0.5  0.5 0.5];
ub=[ub 1    12   20 20  20    3    log(1.1)  3   2    2   2];

if contains(modelName, '_noIns1')
    lb([28 31])=[];
    ub([28 31])=[];
elseif contains(modelName, '_noIns2')
    lb([29 32])=[];
    ub([29 32])=[];
elseif  contains(modelName, '_noIns3')
    lb([27 30])=[];
    ub([27 30])=[];
end
nParams=length(lb);

%% setup stimulus
stimulus=table();
ins=unique([expData.FFA.Ins; expData.Glycerol.Ins]);
stimulus.Ins=[ins; 0];
stimulus.Iso=[0.01*ones(height(stimulus)-1,1); 0]; % 10nm = 0.01 µM

%% Setup optimization
opts.ndiverse     = 500; %'auto'; %100; %500; %5; %
opts.maxtime      = 750; % In cess this option will be overwritten
opts.maxeval      = 1e7;
opts.log_var      = [];

opts.local.solver = 'dhc'; %'dhc'; %'fmincon'; %'nl2sol'; %'mix'; %
opts.local.finish = opts.local.solver;
opts.local.bestx = 0;
opts.local.tol = 2;
opts.local.iterprint = 1;

opts.dim_refset   = 'auto'; %

if(strcmp(opts.local.solver,'fmincon'))
    opts.local.use_gradient_for_finish = 1; %DW: provide gradient to fmincon
else
    opts.local.use_gradient_for_finish = 0; %DW: provide gradient to fmincon
end
opts.local.check_gradient_for_finish = 0; %DW: gradient checker

end

