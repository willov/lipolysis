function [] = PlotInVivo(data, sim, id, titleStr, vars, posOverride)
if nargin<4; titleStr=''; end
if nargin <5 || isempty(vars)
  vars=sim.Properties.VariableNames;
  vars(ismember(vars,'Time'))=[];
end

set(figure(id), 'outerposition',[0 0 2560 1440], 'PaperType','a4')

nPlots=round(length(vars)/2);

if nargin>5 && ~isempty(posOverride)
  m=3; n=2; pos=posOverride;
else
  if contains(titleStr, 'Glycerol')
    m=3; n=2; pos=[1 3 5];
  else
    m=2; n=2; pos=[1 2];
  end
end
for i=1:length(vars)
  switch vars{i}
    case 'Fig1', c=[86 180 233]/256; plt=1; leg='adr';
    case 'Fig1Diab', c=[213 94 0]/256; plt=1; leg='adr (T2D)';
    case 'Fig2EpiDiab', c=[213 94 0]/256; plt=pos(2); leg='adr (T2D)';
    case 'Fig2Epi',  c = [86 180 233]/256;  plt=pos(2); leg='adr';
    case 'Fig2Iso',  c = [170 50 177]/256; plt=3; leg='iso';
    case 'Fig3Epi', c = [86 180 233]/256; plt=pos(end); leg='adr';
    case 'Fig3EpiPhe', c = [0 103 28]/256; plt=pos(end); leg='adr + phe';
  end
  if isfield(data,'InVivo')
    PlotExperiment(m,n,plt,data.InVivo.(vars{i}),sim(:,{'Time',vars{i}}), leg, c)
  else
    PlotExperiment(m,n,plt,[],sim(:,{'Time',vars{i}}), leg, c)
  end
end

for i = 1:nPlots
  subplot(m,n,pos(i))
  PlotInput(m,n,pos(i), titleStr)
end
end


function [] = PlotExperiment(m,n,plt, data, sim, leg, c)
subplot(m,n,plt)
hold on
x=sim{:,1}-15;
y=sim{:,2};

if size(y,2)>1 && isequal(y(:,2), y(:,3))
  y=y(:,2);
end

if size(y,2)>1
  h = fill([x; flipud(x)], [y(:,2); flipud(y(:,3))],c,'HandleVisibility','off');
  set(h,'facealpha',0.5,'edgealpha',0);
end

plot(x, y(:,1), 'color', c,'linewidth', 2.5);
hold on
if  isfield(data,'Mean')
  errorbar(data.Time-15, data.Mean, data.SEM, 'o','MarkerFaceColor', 'auto', 'linewidth', 2.5,'capsize',12, 'color', c, 'HandleVisibility','off')
end
lgd=legend();
lgd.String{end}=leg;
end


function [] = PlotInput(m,n,plt, titleStr)
subplot(m,n,plt)
lgd=legend();
lgd0=lgd.String;

axsPatch=findall(gca,'type','patch');
axsLine=findobj(gca,'type','line');
axsErrorbar = findall(gca, 'type', 'errorbar');
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

plot([15 45]-10,[minimum-0.06*(maximum-minimum)  minimum-0.06*(maximum-minimum)],'color',[0.85 0.85 0.85],'linewidth',3)
plot([45 75]-10,[minimum-0.03*(maximum-minimum)  minimum-0.03*(maximum-minimum)],'color',[0.7 0.7 0.7],'linewidth',3)

h=   plot([75 370]-10,[minimum-0.1*(maximum-minimum)  minimum-0.1*(maximum-minimum)],'color',[0 0 0],'linewidth',3);
if plt==1
  set(h,'LineStyle','none')
end
plot([15 45]+180-10,[minimum-0.06*(maximum-minimum)  minimum-0.06*(maximum-minimum)],'color',[0.85 0.85 0.85],'linewidth',3)
plot([45 75]+180-10,[minimum-0.03*(maximum-minimum)  minimum-0.03*(maximum-minimum)],'color',[0.7 0.7 0.7],'linewidth',3)


axis('tight');
a=axis;
axis([a(1), a(2), minimum-0.175*(maximum-minimum), a(4)])
xlabel('Time (min)')
if strcmp(titleStr,'FFA')
  ylabel({'Released FA','% of max response'})
else
ylabel(sprintf('Released %s (\\muM)', titleStr))
end
set(gca,'FontSize', 20)
lgd.String=lgd0;

end







