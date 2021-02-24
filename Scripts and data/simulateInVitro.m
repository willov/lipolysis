function [responses, steady,  responsesRaw, basal] = simulateInVitro(model, params, expInd, diab,  stimulus, plotTS)
if nargin<6 || isempty(plotTS); plotTS=false; end

params(1)=0; % The parameter for drift in in vivo data
params(expInd)=0;% The parameter for clearance of glycerol/FFA

[response, ~, steady, sims]=simulate_dose_response(model,params ,0,stimulus);
responsesRaw.Normal=response;
vars={'Glycerol','FFA','HSL','PKB'};
basal=response(1,vars);
response{:,vars}=response{:,vars}./basal{1,vars}*100;

responses.Normal=response;
if plotTS
    figure(20)
    c=[jet(size(sims,2)-1); 0 0 0];
    variables=IQMvariables(model);
    observables=strrep(variables(contains(variables,'y_')),'y_',''); % Find variables containing "y_", and then remove that string from those variables.
    
    for i=1:size(sims,2)
        for j=1:length(observables)
            subplot(3,3,j)
            hold on
            plot(sims(i).time, sims(i).variablevalues(:,j),'color',c(i,:))
            title(observables{j})
        end
    end
end

if diab~=0
    responses.Normal.Reesterification = 100 * (3*responses.Normal.Glycerol - responses.Normal.FFA)./(3*responses.Normal.Glycerol);
    [response]=simulate_dose_response(model,params, diab,stimulus);
    responsesRaw.Diabetes=response;
    response{:,vars}=response{:,vars}./basal{1,vars}*100;
    responses.Diabetes=response;
    responses.Diabetes.Reesterification = 100 * (3*responses.Diabetes.Glycerol - responses.Diabetes.FFA)./(3*responses.Diabetes.Glycerol);
    
end

end
function [response, basal, steady, sims] = simulate_dose_response(model, params, diab,stimulus)

variables=IQMvariables(model);
observables=strrep(variables(contains(variables,'y_')),'y_',''); % Find variables containing "y_", and then remove that string from those variables.
inputs=variables(~contains(variables,'y_'));
response=array2table(nan(height(stimulus),length(observables)),'VariableNames',observables);
IR=array2table(nan(height(stimulus),length(inputs)),'VariableNames',inputs);
washInd=ismember(IQMstates(model),{'Gly','FFA'});

sim0=model([0 30000-10 30000],[],[params  diab 0 0 0 0]); % simulate rest before start of experiment
steady = abs(sim0.statevalues(end-1,1:end-2)-sim0.statevalues(end,1:end-2)); % ignores FFA and Glycerol released since they are washed away
steady = sum(steady(steady>1e-6))*1e9;
ic=sim0.statevalues(end,:);
ic(washInd)=0; % wash away extracellular stuff

for i = fliplr(1:height(stimulus))
    
    sim0=model((0:15)',ic,[params  diab 0 0 0 0]); %simulate 15 minutes without stimulus
    
    simIncub=model((0:30)',sim0.statevalues(end,:),[params diab 0 0 0 0]);% Others
    simIns=model((simIncub.time(end):simIncub.time(end)+15)',simIncub.statevalues(end,:),[params diab 0 0 0 stimulus.Ins(i)]);%
    simIso=model((simIns.time(end):simIns.time(end)+10)',simIns.statevalues(end,:),[params diab 0 0 stimulus.Iso(i) stimulus.Ins(i)]);
    
    response{i,:}=simIso.variablevalues(end,1:length(observables));
    
    IR{i,:}=simIso.variablevalues(end,end-length(inputs)+1:end);
    sims(i)=StructFieldCat(simIncub, [simIns simIso]);
    if i==1
        basal=simIso.statevalues(end,:);
    end
end
response=[stimulus(:,1:2) IR response];
end


% end