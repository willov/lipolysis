function [responses, scale, responsesRaw]=SimulateInVivo(params, model, data, diab)
ins=data.Fig1.Ins;

if nargin <4, diab= 1; end %not used in most simulations. 1 = no diabition

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

try
    init = model([0 20000-10 20000],[],[params diab 0 0 0 ins],simOptions);
catch err
    cost = 1e90;
    responsesRaw=[];
    responses=[];
    return
end
ic = init.statevalues(end,:);
responsesRaw.Gly=table();
responsesRaw.Gly.Time=(data.Fig1.Time(1):data.Fig1.Time(end))';
responsesRaw.FFA=table();
responsesRaw.FFA.Time=(data.Fig1.Time(1):data.Fig1.Time(end))';

ind = ismember(data.Fig1.Time(1):data.Fig1.Time(end), data.Fig1.Time);

Fig1 = SimulateExperiment(params, diab, model, data.Fig1, ic, [1e3 0 0]);
yh=Fig1(ind, 1);
y=data.Fig1.Mean;
SEM=data.Fig1.SEM;
scale1=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig1=Fig1(:,1);
responsesRaw.FFA.Fig1=Fig1(:,2);

Fig2Epi = SimulateExperiment(params, diab, model, data.Fig2Epi, ic, [1e3 0 0]);
Fig2Iso = SimulateExperiment(params, diab, model, data.Fig2Iso, ic, [0 1e2 0]);
yh=[Fig2Epi(ind, 1); Fig2Iso(ind, 1)];
y=[data.Fig2Epi.Mean; data.Fig2Iso.Mean];
SEM=[data.Fig2Epi.SEM; data.Fig2Iso.SEM];
scale2=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig2Epi=Fig2Epi(:,1);
responsesRaw.Gly.Fig2Iso=Fig2Iso(:,1);
responsesRaw.FFA.Fig2Epi=Fig2Epi(:,2);
responsesRaw.FFA.Fig2Iso=Fig2Iso(:,2);

Fig3Epi = SimulateExperiment(params, diab, model, data.Fig3Epi, ic, [1e3 0 0]);
Fig3EpiPhe = SimulateExperiment(params, diab, model, data.Fig3EpiPhe, ic, [1e3 0 1]);
yh=[Fig3Epi(ind, 1); Fig3EpiPhe(ind, 1)];
y=[data.Fig3Epi.Mean; data.Fig3EpiPhe.Mean];
SEM=[data.Fig3Epi.SEM; data.Fig3EpiPhe.SEM];
scale3=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig3Epi=Fig3Epi(:,1);
responsesRaw.Gly.Fig3EpiPhe=Fig3EpiPhe(:,1);
responsesRaw.FFA.Fig3Epi=Fig3Epi(:,2);
responsesRaw.FFA.Fig3EpiPhe=Fig3EpiPhe(:,2);

scale=[scale1 scale2 scale2 scale3 scale3];

responses=responsesRaw;
responses.Gly{:,2:end}=responses.Gly{:,2:end}.*scale;
responses.FFA{:,2:end}=responses.FFA{:,2:end}.*scale;

end

function [response] = SimulateExperiment(params, diab, model,  data, ic, stimulus)
simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

Adr = stimulus(1);
Iso = stimulus(2);
Phe = stimulus(3);
if diab~=1 %increase basal levels of insulin in diabetes
   data.Ins(data.Ins==min(data.Ins)) = min(data.Ins)*2;
end
Ins = min(data.Ins);

t = data.Time;
ind=contains(IQMvariables(model),{'Gly','FFA'});

try
    init = model([0 2000-10 2000],ic,[params diab 0 0 0 Ins],simOptions);
    
    ic = init.statevalues(end,:);
    
    simBasal = model(t(1):t(2)+5,ic,[params diab 0 0 0 Ins],simOptions); %Steady state
    simLow = model(t(2)+5:t(5)+5,simBasal.statevalues(end,:),[params diab Phe Adr Iso Ins],simOptions); % low dose
    simHigh = model(t(5)+5:t(8)+5,simLow.statevalues(end,:),[params diab Phe Adr*10 Iso*10 Ins],simOptions); % high dose
    
    Ins = max(data.Ins); %start infusion of insulin
    simIns = model(t(8)+5:t(16)+5,simHigh.statevalues(end,:),[params diab 0 0 0 Ins],simOptions); % insulin only
    simInsLow = model(t(16)+5:t(19)+5,simIns.statevalues(end,:),[params diab Phe Adr Iso Ins],simOptions); % low dose
    simInsHigh = model(t(19)+5:t(22)+5,simInsLow.statevalues(end,:),[params diab Phe Adr*10 Iso*10 Ins],simOptions); %high dose
    simIns2 = model(t(22)+5:t(30),simInsHigh.statevalues(end,:),[params diab 0 0 0 Ins],simOptions); %insulin only

    response = [simBasal.variablevalues(:,ind);
        simLow.variablevalues(2:end,ind);
        simHigh.variablevalues(2:end,ind);
        simIns.variablevalues(2:end,ind);
        simInsLow.variablevalues(2:end,ind);
        simInsHigh.variablevalues(2:end,ind);
        simIns2.variablevalues(2:end,ind)];
    
    
catch err
    response= repmat(1e90,fliplr(size(t(1):t(end))));
    return
end



end