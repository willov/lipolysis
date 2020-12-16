function [] = PlotInVitro(responses, expData, experiments, id, xVar, posOverride)
%PLOT_DOSE_RESPONSE Summary of this function goes here
%   Detailed explanation goes here

if nargin<5 || isempty(xVar)
    xVar='Ins';
end
if nargin<3 || isempty(experiments)
  experiments = fieldnames(responses);
end

if nargin > 5 && ~isempty(posOverride)
  pos = posOverride;
  m=3;
  n=2;
else
  if width(responses.Normal)==4
    m=2; n=2; pos=[1 2];
  elseif width(responses.Normal)==3
    m=2; n=2; pos=[1];
  else
    m=3; n=2; pos=[2 4 6];
  end
end
if nargin > 5 && ~isempty(posOverride)
  pos = posOverride;
end
figure(id);
for i =1:length(experiments)
    switch experiments{i}
        case 'Normal'
            c=[2 64 167]/256;
        case 'Wort'
            c=[0.2344    0.6484    0.2891];
        case 'Torin'
            c=[129, 25, 133]/256;
        case 'Akti'
            c=[168, 88, 42]/256;
        case 'Diabetes'
            c=[230 159 0]/256;
        otherwise
            c=[0 0 0];        
    end
    
    if strcmp(experiments{i}, 'Normal')
        dontPlotVars=[];
    elseif strcmp(experiments{i}, 'Diabetes')
        response=responses.(experiments{i});
        dontPlotVars={};
    else
        dontPlotVars={'IR','IR_2','IR_3'};
    end
    
    response=responses.(experiments{i});
    if ~isempty(dontPlotVars)
        response{:,dontPlotVars}=nan(size(response{:,dontPlotVars}));
    end
    
    if isfield(expData,experiments{i})
        data=expData.(experiments{i});
    else
        data=[];
    end
    PlotExperiment(m,n,pos,response, data, c, xVar)
    
end
PlotInput(m,n,pos,responses.Normal)
set(figure(id), 'outerposition',[0 0 2560 1440], 'PaperType','a4')
end

function []=PlotExperiment(m,n,pos,response, data, c, xVar)

set(0,'DefaultLineLineWidth',2)

for j=3:width(response)
    variable=response.Properties.VariableNames{j};
    if isfield(data,variable)
        PlotSubplot(data.(variable), response(:,[1 2 j]),m,n,pos(j-2), c, xVar)
    else
        PlotSubplot([], response(:,[1 2 j]),m,n,pos(j-2), c, xVar)
    end
end
end

function []=PlotSubplot(data, sim, m,n,ind,c, xVar)
x=sim.(xVar)*1e-9;
if x(1)==0
    x(1)=x(2)/4;
end
if x(end)==0
    x(end)=x(end-1)*10;
end

subplot(m,n,ind)

y=sim{:,end};
if size(y,2)>1 && isequal(y(:,2), y(:,3))
    y=y(:,2);
end

if size(y,2)==3
    h  = fill([x(1)/1.5 x(1)/1.5 x(1)*1.5 x(1)*1.5],[y(1, 2) y(1, 3) y(1, 3) y(1, 2)  ], c);
    set(h,'facealpha',0.5,'edgealpha',0);
    hold on
    h2 = fill([x(2:end-1); flipud(x(2:end-1))], [y(2:end-1,2); flipud(y(2:end-1,3))],c);
    set(h2,'facealpha',0.5,'edgealpha',0);
    h3  = fill([x(end)/1.5 x(end)/1.5 x(end)*1.5 x(end)*1.5],[y(end, 2) y(end, 3) y(end, 3) y(end, 2)  ], c);
    set(h3,'facealpha',0.5,'edgealpha',0);
elseif size(y,2)>3
    plot([x(1)/1.5 x(1)*1.5],[y(1, :); y(1, :)],'-','Color',c,'linewidth', 2.5);
    hold on
    semilogx(x(2:end-1),y(2:end-1,:),'-','Color',c)
    plot([x(end)/1.5 x(end)*1.5],[y(end, :); y(end, :)],'-','Color',c,'linewidth', 2.5);
end
if any(ismember(size(y,2),[1 3]))
    plot([x(1)/1.5 x(1)*1.5],[y(1, 1) y(1, 1)],'-','Color',c,'linewidth', 2.5);
    hold on
    
    plot(x(2:end-1), y(2:end-1,1), 'color', c,'linewidth', 2.5);
    plot([x(end)/1.5 x(end)*1.5],[y(end, 1) y(end, 1)],'-','Color',c,'linewidth', 2.5);
end

if ~isempty(data)
    insData=data.Ins*1e-9;
    insData([1 end])=[x(1) x(end)]; % adjust for zeros in logscale (use same as in simulation)
    errorbar(insData, data.Mean, data.SEM,'o','MarkerFaceColor', 'auto', 'linewidth', 2.5,'capsize',12,'Color',c)
end
end

function [] = PlotInput(m,n,pos, response)

for i = 1:width(response)-2
    subplot(m,n,pos(i))
    axis('tight');
    
    xlabel('[Insulin] (M)')
    if contains(response.Properties.VariableNames{i+2},{'HSL', 'PKB'})
        ylabel({sprintf('Phosphorylation of %s,',response.Properties.VariableNames{i+2}); '% of iso only'})
    else
        ylabel({sprintf('Released %s,',response.Properties.VariableNames{i+2}); '% of iso only'})
    end
    
    axsPatch=findobj(gca,'type','patch');
    axsLine=findobj(gca,'type','line');
    axsErrorbar = findobj(gca, 'type', 'errorbar');
    YData=[];
    if ~isempty(axsPatch)
        YData=[YData horzcat(axsLine.YData)];
    end
    if ~isempty(axsLine)
        YData=[YData vertcat(axsPatch.YData)'];
    end
    if ~isempty(axsErrorbar)
        YData=[YData horzcat(axsErrorbar.YData)-horzcat(axsErrorbar.YNegativeDelta) horzcat(axsErrorbar.YData)+horzcat(axsErrorbar.YPositiveDelta)];
    end
    minimum = min(YData);
    maximum = max(YData);
    
    x=response.Ins*1e-9;
    x(1)=x(2)/4/1.5;
    plot(x([1 end-1]), [minimum-0.03*(maximum-minimum) minimum-0.03*(maximum-minimum)],'color', [0.7 0.7 0.7],'linewidth', 3) % iso stimulus bar
    plot(x([2 end-1]), [minimum-0.1*(maximum-minimum) minimum-0.1*(maximum-minimum) ],'color', [0 0 0],'linewidth', 3) % ins stimulus bar
    set(gca,'XScale','log', 'XMinorTick','on','FontSize', 20)
    
    axis('tight');
    a=axis;
    axis([a(1), a(2), minimum-0.175*(maximum-minimum), a(4)])
    box off
end
end