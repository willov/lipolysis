function [] = PlotPL(modelName, useFig2EpiData)


if nargin<1; modelName='lipolysis'; end
if nargin<2; useFig2EpiData=1; end

    if useFig2EpiData
        baseFolder='./Parameter sets';
    else
        baseFolder='./Parameter sets (with alternative epi data)';
    end
folder= 'Parameter sets/PL/';

[model, data, lb, ub, nParams, expInd, stimulus, dgf] = Init(modelName, 0);
pNames = IQMparameters(model);
rejectlb = lb;
rejectub = ub;
for i = 1:nParams
    rejectlb(i) = NonIdentifiabilityBound(pNames{i}, lb(i), 'min');
    rejectub(i) = NonIdentifiabilityBound(pNames{i}, ub(i), 'max');
end

files=dir(sprintf('%s/**/*%s*.mat', folder, modelName));
optParams=[];
for i = fliplr(1:length(files))
    load([files(i).folder '/' files(i).name],'optParam');
    optParams(i,:)=optParam;
end
optParams=unique(optParams,'rows');
%% Print bounds
maxP = max(optParams)';
maxP=[exp(maxP(1:expInd)); maxP(expInd+1:end)];
maxP(27:29) = 10.^maxP(27:29);
minP = min(optParams)';
minP=[exp(minP(1:expInd)); minP(expInd+1:end)];
minP(27:29) = 10.^minP(27:29);
lowerThreshold =[exp(rejectlb(1:expInd)) rejectlb(expInd+1:end)]';
upperThreshold =[exp(rejectub(1:expInd)) rejectub(expInd+1:end)]';

params=table(pNames(1:32),lowerThreshold, upperThreshold, minP,maxP,'VariableNames',{'Parameters', 'Lower_threshold', 'Upper_threshold','Min','Max'});
disp('Boundry for parameter uncertainty, presented in Table S1.')
disp(params)
%%
x = 2*(optParams-rejectlb) ./ (rejectub-rejectlb)-1;

minx = min(x);
maxx=max(x);

figure(2)
clf
plot([0 nParams],[1 1],'color',[0.6350 0.0780 0.1840])
hold on
plot([0 nParams],[-1 -1],'color', [0.6350 0.0780 0.1840])

%%
center = (maxx+minx)/2;
bound = (maxx-minx)/2;
h = errorbar(1:size(x,2),center, bound,'.','linewidth',3, 'color', 'k');
h.CapSize = 20;
xticks(1:size(x,2))
xticklabels(IQMparameters(model))
set(gca,'TickLabelInterpreter','none', 'fontsize',30)
xtickangle(45)
ylim([-1.05 1.05])
box off
ylabel({'Non-rejected parameter values,', 'relative to the threshold for non-identifiability'})
hold on

%% Plot the bounds relative to the non-identifiability threshold
relativelb = 2*(lb-rejectlb) ./ (rejectub-rejectlb)-1;
relativeub = 2*(ub-rejectlb) ./ (rejectub-rejectlb)-1;
for i = 1:length(relativelb)
    plot([i-0.3 i+0.3], [relativelb(i) relativelb(i)], 'color', [0.9290 0.6940 0.1250])
    plot([i-0.3 i+0.3], [relativeub(i) relativeub(i)], 'color', [0.9290 0.6940 0.1250])
end

%% Plot the best solution
load(FindBestParametersFile(baseFolder, 1, [modelName ', opt-eSS']), 'optParam')
x = 2*(optParam-rejectlb) ./ (rejectub-rejectlb)-1;
plot(x,'x','MarkerSize',15, 'color',[0 0.4470 0.7410]);

%% Set size and save figure
set(gcf, 'outerposition',[0 0 2560 1440], 'PaperType','a4')


end
