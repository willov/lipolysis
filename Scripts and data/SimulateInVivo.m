function [responses, penalty,  scale, responsesRaw]=SimulateInVivo(params, model, data, diab, tol)
ins=data.Fig1.Ins;

if nargin<4, diab= 0; end %not used in most simulations. 0 = no diabition
if nargin<5, tol =0;end

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

try
    init = model([0 30000-10 30000],[],[params diab 0 0e-4 0 ins],simOptions);
    
    steady = abs(init.statevalues(end-1,:)-init.statevalues(end,:));
    steady = sum(steady(steady>1e-6))*1e9;
    
catch err
    cost = 1e90;
    responsesRaw=[];
    responses=[];
    steady=1e20;
    return
end

ic = init.statevalues(end,:);
responsesRaw.Gly=table();
responsesRaw.Gly.Time=(data.Fig1.Time(1):data.Fig1.Time(end))';
responsesRaw.FFA=table();
responsesRaw.FFA.Time=(data.Fig1.Time(1):data.Fig1.Time(end))';

ind = ismember(data.Fig1.Time(1):data.Fig1.Time(end), data.Fig1.Time);

[Epi, simEpi] = SimulateExperiment(params, diab, model, data.Fig1, ic, [1 0 0]);
yh=Epi(ind, 1);
yh = [ones(size(yh)) yh];
y=data.Fig1.Mean;
SEM=data.Fig1.SEM;
scale1=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig1=Epi(:,1);
responsesRaw.FFA.Fig1=Epi(:,2);
penaltyEpi = Penalty(simEpi.statevalues, tol);

[EpiIns, simEpiIns] = SimulateExperiment(params, diab, model, data.Fig2Epi, ic, [1 0 0]);
[IsoIns, simIsoIns] = SimulateExperiment(params, diab, model, data.Fig2Iso, ic, [0e-4 0.1 0]);
yh=[EpiIns(ind, 1); IsoIns(ind, 1)];
yh = [ones(size(yh)) yh];
y=[data.Fig2Epi.Mean; data.Fig2Iso.Mean];
SEM=[data.Fig2Epi.SEM; data.Fig2Iso.SEM];
scale2=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig2Epi=EpiIns(:,1);
responsesRaw.Gly.Fig2Iso=IsoIns(:,1);
responsesRaw.FFA.Fig2Epi=EpiIns(:,2);
responsesRaw.FFA.Fig2Iso=IsoIns(:,2);
penaltyEpiIns = Penalty(simEpiIns.statevalues, tol); %3-4 is equal to the alpha receptor, which does not change if no epi is present
penaltyIsoIns = Penalty(simIsoIns.statevalues(:,[1:2 5:end]), tol);

[EpiInsPhe, simEpiPhe] = SimulateExperiment(params, diab, model, data.Fig3EpiPhe, ic, [1 0 1]);
yh=[EpiIns(ind, 1); EpiInsPhe(ind, 1)];
yh = [ones(size(yh)) yh];
y=[data.Fig3Epi.Mean; data.Fig3EpiPhe.Mean];
SEM=[data.Fig3Epi.SEM; data.Fig3EpiPhe.SEM];
scale3=lscov(yh,y, SEM.^-2);
responsesRaw.Gly.Fig3Epi=EpiIns(:,1);
responsesRaw.Gly.Fig3EpiPhe=EpiInsPhe(:,1);
responsesRaw.FFA.Fig3Epi=EpiIns(:,2);
responsesRaw.FFA.Fig3EpiPhe=EpiInsPhe(:,2);
penaltyEpiInsPhe = Penalty(simEpiPhe.statevalues, tol);
scale=[scale1 scale2 scale2 scale3 scale3];

responses=responsesRaw;
responses.Gly{:,2:end}=responses.Gly{:,2:end}.*scale(2,:) + scale(1,:);
responses.FFA{:,2:end}=responses.FFA{:,2:end}.*scale(2,1) + scale(1,1);

penalty = steady+penaltyEpi+penaltyEpiIns+penaltyIsoIns+penaltyEpiInsPhe; % Force the simulations to have reasonable (unscaled) values.
end

function [penalty] = Penalty(sim, tol)
if nargin <2, tol = 0; end
if any(max(sim)-min(sim) <= 0.1-tol)
    penalties = 1./(max(sim)-min(sim));
    penalty = sum(penalties(isfinite(penalties)))+100;
else
    penalty = 0;
end
end

function [response, sim] = SimulateExperiment(params, diab, model,  data, ic, stimulus)
simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

Adr = stimulus(1);
Iso = stimulus(2);
Phe = stimulus(3);
if diab~=0 %increase basal levels of insulin in diabetes
    data.Ins(data.Ins==min(data.Ins)) = min(data.Ins)*2;
end
Ins = min(data.Ins);

t = data.Time;
ind=contains(IQMvariables(model),{'Gly','FFA'});

try
    simBasal = model(t(1):t(2)+5,ic,[params diab 0 0e-4 0 Ins],simOptions); %Steady state
    simLow = model(t(2)+5:t(5)+5,simBasal.statevalues(end,:),[params diab Phe Adr Iso Ins],simOptions); % low dose
    simHigh = model(t(5)+5:t(8)+5,simLow.statevalues(end,:),[params diab Phe Adr*10 Iso*10 Ins],simOptions); % high dose
    
    Ins = max(data.Ins); %start infusion of insulin
    simIns = model(t(8)+5:t(16)+5,simHigh.statevalues(end,:),[params diab 0 0e-4 0 Ins],simOptions); % insulin only
    simInsLow = model(t(16)+5:t(19)+5,simIns.statevalues(end,:),[params diab Phe Adr Iso Ins],simOptions); % low dose
    simInsHigh = model(t(19)+5:t(22)+5,simInsLow.statevalues(end,:),[params diab Phe Adr*10 Iso*10 Ins],simOptions); %high dose
    simIns2 = model(t(22)+5:t(30),simInsHigh.statevalues(end,:),[params diab 0 0e-4 0 Ins],simOptions); %insulin only
    
    response = [simBasal.variablevalues(:,ind);
        simLow.variablevalues(2:end,ind);
        simHigh.variablevalues(2:end,ind);
        simIns.variablevalues(2:end,ind);
        simInsLow.variablevalues(2:end,ind);
        simInsHigh.variablevalues(2:end,ind);
        simIns2.variablevalues(2:end,ind)];
    sim = StructsFieldCat([simBasal;simLow;simHigh;simIns;simInsLow;simInsHigh;simIns2]); %combines multiple structs to one.
    
catch err
    response= repmat(1e90,fliplr(size(t(1):0.01:t(end))));
    return
end



end