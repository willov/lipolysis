function [responses, responsesRaw, basal] = simulateDR(model, params, expInd, doDiabetes, diab,  stimulus, plotTS)
if nargin<7 || isempty(plotTS); plotTS=false; end

params(1)=0; % The parameter for drift in in vivo data
params(expInd)=0;% The parameter for clearance of glycerol/FFA

[response, ~, sims]=simulate_dose_response(model,params ,[1],stimulus);
responsesRaw.Normal=response;
vars={'cAMP','Glycerol','FFA','HSL','PKB'};
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

if doDiabetes %Only do the diabetes simulation, if the diabets parameters is smaller than one, with 1e-1 tolerance due to precision
    [response]=simulate_dose_response(model,params, diab,stimulus);
    response{:,vars}=response{:,vars}./basal{1,vars}*100;
    responses.Diabetes=response;
end

end
function [response, basal, sims] = simulate_dose_response(model, params, diab,stimulus)
% Steady state
variables=IQMvariables(model);
observables=strrep(variables(contains(variables,'y_')),'y_',''); % Find variables containing "y_", and then remove that string from those variables.
        inputs=variables(~contains(variables,'y_'));
        response=array2table(nan(height(stimulus),length(observables)),'VariableNames',observables);
        IR=array2table(nan(height(stimulus),length(inputs)),'VariableNames',inputs);
        washInd=ismember(IQMstates(model),{'Gly','FFA'});
        % sims(height(stimulus))=struct();
        for i = fliplr(1:height(stimulus))
            sim0=model([0 10000],[],[params  diab 0 0 0 0]); % simulate rest before start of experiment
            ic=sim0.statevalues(end,:);
            ic(washInd)=0; % wash away extracellular stuff
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