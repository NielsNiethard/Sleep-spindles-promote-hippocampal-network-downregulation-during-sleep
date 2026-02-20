% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
%% Triplet Analysis
close all
clear

addpath('/path/to/software/Mass_Univariate_ERP_Toolbox-master')

Bsl = [2; 3];
ClusterName{1,1} = 'SO';
ClusterName{2,1} = 'Spindle';
ClusterName{3,1} = 'SO+Spindle';
ClusterName{4,1} = 'SO+Spindle2';
PlotColors{1,1} = [1 0 0];
PlotColors{2,1} = [0 1 0];
PlotColors{3,1} = [0 0 1];
PlotColors{4,1} = [0 0 1];
MicroaeousalMerge =1;
DirData = '/path/to/Results/CaData/';
MinEpisodeDuration = 40;
Files = dir(DirData);

for iFile = 3: size(Files,1)
    FileName{iFile-2,1} = Files(iFile,1).name;
end

for iCluster = 1:4
    switch iCluster
        case 1
            DoCoupledEvents = 0;
            DoSpindleEvents = 0;
        case 2
            DoCoupledEvents = 0;
            DoSpindleEvents = 1;
        case 3
            DoCoupledEvents = 1;
            DoSpindleEvents = 1;
        case 4
            DoCoupledEvents = 1;
            DoSpindleEvents = 0;
    end
    
    GAAllDS = [];
    GAAll = [];
    ActiveCellsFractionAll = [];
    EventDurationAll = [];
    CaDuringEvents =[];
    
    for iFile = 1:length(FileName)
        load(strcat(DirData,FileName{iFile,1}),'Events','CaDataAll','SlstAll','FrameRate');
        SpiSO = nan(1,size(Events.Slo,2));
        
        for iEvent = 1: size(Events.Slo,2)
            if ~isempty(find(Events.Slo(2,iEvent)<Events.Spi(1,:) & (Events.Slo(2,iEvent) + 1.5) > Events.Spi(1,:)))
                SpiSO(:,iEvent) = [min(find(Events.Slo(2,iEvent)<Events.Spi(1,:) & (Events.Slo(2,iEvent) + 1.5) > Events.Spi(1,:)))];
            end
        end
        SpiSO(isnan(SpiSO))=[];
        SpiSO = unique(SpiSO);
        
        SOSpi = nan(1,size(Events.Spi,2));
        for iEvent = 1: size(Events.Spi,2)
            if ~isempty(find(Events.Spi(2,iEvent)>Events.Slo(2,:) & Events.Spi(2,iEvent) < (Events.Slo(2,:)+1.5)))
                SOSpi(:,iEvent) = [max(find(Events.Spi(2,iEvent)>Events.Slo(2,:) & Events.Spi(2,iEvent) < (Events.Slo(2,:)+1.5)))];
            end
        end
        SOSpi(isnan(SOSpi))=[];
        SOSpi = unique(SOSpi);
        
        if DoCoupledEvents == 1
            if DoSpindleEvents ==  1
                CurrentEvents = Events.Spi(:,SpiSO);
            else
                CurrentEvents = [Events.Slo(4,SOSpi);Events.Slo(3,SOSpi)];
            end
        else
            if DoSpindleEvents ==  1
                CurrentEvents = Events.Spi;
                CurrentEvents(:,SpiSO)=[];
            else
                CurrentEvents = [Events.Slo(4,:);Events.Slo(3,:)];
                CurrentEvents(:,SOSpi)=[];
            end
        end
        
        WindowSize = 5;
        TMPCaDataEventGADS = nan(size(CurrentEvents,2),round(WindowSize*31)+round(WindowSize*31)+1);
        TMPCaDataEventGA = nan(size(CurrentEvents,2),round(WindowSize*31)+round(WindowSize*31)+1);
        EventDuration = nan(size(CurrentEvents,2),1);
        TmpCaDuringEvents = nan(size(CurrentEvents,2),1);
        ActiveCellsFraction = nan(size(CurrentEvents,2),1);
        
        for iEvent = 1:size(CurrentEvents,2)
            if round(CurrentEvents(1,iEvent)*31)+round(WindowSize*31) < size(CaDataAll,2) &&...
                    (round(CurrentEvents(1,iEvent)*31)-round(WindowSize*31)) > 0
                
                TMPAct =nansum(~isnan(CaDataAll(:,round(CurrentEvents(1,iEvent)*31)-round(WindowSize*31):...
                    round(CurrentEvents(1,iEvent)*31)+round(WindowSize*31))),1); % active cells per timepoint
                
                TMPAmp =nansum((CaDataAll(:,round(CurrentEvents(1,iEvent)*31)-round(WindowSize*31):...
                      round(CurrentEvents(1,iEvent)*31)+round(WindowSize*31)))); %sum amplitude per timepoint
                      
                TMPCaDataEventGA((iEvent),:) = TMPAmp./TMPAct;
                
                EventDuration(iEvent,1) = (CurrentEvents(2,iEvent)-CurrentEvents(1,iEvent));
                
                ActiveCellsFraction(iEvent,1) = 100 * sum(~isnan(nanmean(CaDataAll(:,round(CurrentEvents(1,iEvent)*31):...
                    round(CurrentEvents(2,iEvent)*31)),2)))/size(CaDataAll,1);
            end
        end
        
        ActiveCellsFraction(isnan(ActiveCellsFraction)) = 0;
        EventDuration(isnan(nanmean(TMPCaDataEventGA,2)),:) = [];
        TMPCaDataEventGA(isnan(nanmean(TMPCaDataEventGA,2)),:) = [];
        
        GAAll = [GAAll; TMPCaDataEventGA];
        EventDurationAll = [EventDurationAll; EventDuration];
        ActiveCellsFractionAll = [ActiveCellsFractionAll; ActiveCellsFraction];
    end
    
    Output.GAAll{iCluster} =GAAll;
    Output.EventDurationAll{iCluster} =EventDurationAll;
    Output.ActiveCellsFraction{iCluster} =ActiveCellsFractionAll;
end

for iCluster = 1:4
    spikeFreq{iCluster} = Output.GAAll{iCluster};
    spikeFreq{iCluster}(isnan(spikeFreq{iCluster}))=0;
    spikeFreq{iCluster} = movmean(spikeFreq{iCluster},10,2);
end
%%
DSFactor = 1;

Output.CaDSFreq = spikeFreq;

for iCluster = 1:4
    hFig = figure(1000);
    if iCluster == 1
        tiledlayout(4,6)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 1000 700])
        set(gcf,'color','white')
        set(gcf,'name',' Grand average','NumberTitle','off')
    end
    nexttile((iCluster-1)*6+1)
    
    Output.CaDSFreq{iCluster}(isnan(Output.CaDSFreq{iCluster})) = 0;
    
    stdshade(Output.CaDSFreq{iCluster}(:,:)-nanmean(nanmean(Output.CaDSFreq{iCluster}(:,Bsl(1,1)*31/DSFactor:Bsl(2,1)*31/DSFactor),2),1),0.5,PlotColors{iCluster,1});
    hold all
    
    xlim([62/DSFactor 248/DSFactor])
    xline(155/DSFactor)
    ylim([0 0.2])
    Labels = {'-3','-2','-1','0','1','2','3'};
    set(gca, 'XTick', 62/DSFactor:31/DSFactor:248/DSFactor, 'XTickLabel', Labels);
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca, 'box', 'off')
    title(ClusterName{iCluster,1})
    
    PValue = [];
    for iTimepoint = 1:size(Output.CaDSFreq{iCluster},2)
        [h,p] = ttest(Output.CaDSFreq{iCluster}(:,iTimepoint)-nanmean(nanmean(Output.CaDSFreq{iCluster}(:,Bsl(1,1)*31/DSFactor:Bsl(2,1)*31/DSFactor),2),1));
        PValue(iTimepoint) = p;
    end
    
    PValue(PValue > 0.05) = nan;
    PValuePlot = PValue;
    PValuePlot(PValue < 0.05) = -0.0002;
    hold on
    plot(1:size(Output.CaDSFreq{iCluster},2),PValuePlot, LineWidth= 4, Color= PlotColors{iCluster,1});
    PValuePlot(PValue > 0.01) = nan;
    plot(1:size(Output.CaDSFreq{iCluster},2),PValuePlot, LineWidth= 8, Color= PlotColors{iCluster,1});
end

for iCluster = 1:4
    for iEvent = 1: size(Output.EventDurationAll{iCluster},1)
        Output.CaDuringEventDuration{iCluster}(iEvent,1) = nanmean(Output.CaDSFreq{iCluster}(iEvent,round(WindowSize*31/DSFactor):...
            round(WindowSize*31/DSFactor)+round(Output.EventDurationAll{iCluster} * 31/DSFactor)),2);
    
        Output.CaBeforeEvent{iCluster}(iEvent,1) = nanmean(Output.CaDSFreq{iCluster}(iEvent,round(WindowSize*31/DSFactor)+(1*31/DSFactor):...
            round(WindowSize*31/DSFactor)+(3*31/DSFactor)),2);   
        Output.CaBaseline{iCluster}(iEvent,1) = nanmean(Output.CaDSFreq{iCluster}(iEvent,round(Bsl(1,1)*31/DSFactor):round(Bsl(2,1)*31/DSFactor)),2);
    end
end

nexttile(2)
for iEvent = 1:4
    b = boxchart(repmat(iEvent,size(Output.CaDuringEventDuration{iEvent},1),1),Output.CaDuringEventDuration{iEvent},'MarkerStyle','none');
    hold all
    b.BoxWidth = 0.75;
    b.BoxFaceColor = PlotColors{iEvent,1};
    b.BoxFaceAlpha = 0.4;
    b.LineWidth = 1.5;
    b.Notch = 'on';
end

ylabel('Mean activity');
title('During Event');
Labels = {'SO','Spindle','SO+spindle','SO+spindle2'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    Data = Output.CaDuringEventDuration{barIndex};
    sorted_data_y = sort(Data, 'ascend');
    sorted_data_x = [];
    [N,edges] = histcounts(Data(:,1),100);
    for iBin = 1: size(N,2)
        if N(1,iBin) >0
            XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
            if N(1,iBin) == 1
                XtmpLocation = 0;
            end
            sorted_data_x = [sorted_data_x XtmpLocation];
        end
    end
    sorted_data_x = sorted_data_x/10+barIndex;
    StandardError(barIndex) = nanstd(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  nanmean(Data);
end
e = errorbar([1:size(StandardError,2)],tmp,StandardError,'.');
clear tmp
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
plotSignificanceBetweenNonParametric(Output.CaDuringEventDuration)

ax = nexttile(3,[1 3]);
[p,tbl,stats] = anovan([Output.CaDuringEventDuration{1};Output.CaDuringEventDuration{2};Output.CaDuringEventDuration{3}],...
    {[ones(size(Output.CaDuringEventDuration{1},1),1); ones(size(Output.CaDuringEventDuration{2},1),1) *2; ones(size(Output.CaDuringEventDuration{3},1),1) *3]},'display','off');
VariableNames = tbl(1,:);
ColumnNames = tbl(2:end,1);
tbl = cell2table(tbl);
tbl.Properties.VariableNames = VariableNames;
tbl(1,:) = [];
table2array(tbl)
axisTable(ax, cell2mat(table2array(tbl(1,[3 6 7]))),{'d.f.','F','Prob>F'},ColumnNames(1,1))
axis off
title('One way ANOVA')
[p, tbl, stats] = kruskalwallis([Output.CaDuringEventDuration{1};Output.CaDuringEventDuration{2};Output.CaDuringEventDuration{3}],...
    [ones(size(Output.CaDuringEventDuration{1},1),1); ones(size(Output.CaDuringEventDuration{2},1),1) *2; ones(size(Output.CaDuringEventDuration{3},1),1) *3],'off');

nexttile(8)
b = bar([nanmean(Output.CaBeforeEvent{1}),...
    nanmean(Output.CaBeforeEvent{2}),...
    nanmean(Output.CaBeforeEvent{3}),...
    nanmean(Output.CaBeforeEvent{4}),]);
b.FaceColor = 'flat';
b.CData([1],:) = PlotColors{1,1};
b.CData([2],:) = PlotColors{2,1};
b.CData([3],:) = PlotColors{3,1};
b.CData([4],:) = PlotColors{4,1};
ylabel('Mean activity');
title('1s Before Event');
Labels = {'SO','Spindle','SO+spindle','SO+spindle2'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    Data = Output.CaBeforeEvent{barIndex};
    sorted_data_y = sort(Data, 'ascend');
    sorted_data_x = [];
    [N,edges] = histcounts(Data(:,1),100);
    for iBin = 1: size(N,2)
        if N(1,iBin) >0
            XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
            if N(1,iBin) == 1
                XtmpLocation = 0;
            end
            sorted_data_x = [sorted_data_x XtmpLocation];
        end
    end
    sorted_data_x = sorted_data_x/10+barIndex;
    StandardError(barIndex) = nanstd(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  nanmean(Data);
end
e = errorbar([1:size(StandardError,2)],tmp,StandardError,'.');
clear tmp
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
plotSignificanceBetweenNonParametric(Output.CaBeforeEvent)

ax = nexttile(9,[1 3]);
[p,tbl,stats] = anovan([Output.CaBeforeEvent{1};Output.CaBeforeEvent{2};Output.CaBeforeEvent{3}],...
    {[ones(size(Output.CaBeforeEvent{1},1),1); ones(size(Output.CaBeforeEvent{2},1),1) *2; ones(size(Output.CaBeforeEvent{3},1),1) *3]},'display','off');
VariableNames = tbl(1,:);
ColumnNames = tbl(2:end,1);
tbl = cell2table(tbl);
tbl.Properties.VariableNames = VariableNames;
tbl(1,:) = [];
table2array(tbl);
axisTable(ax, cell2mat(table2array(tbl(1,[3 6 7]))),{'d.f.','F','Prob>F'},ColumnNames(1,1))
axis off
title('One way ANOVA')
[p, tbl, stats] = kruskalwallis([Output.CaBeforeEvent{1};Output.CaBeforeEvent{2};Output.CaBeforeEvent{3}],...
    [ones(size(Output.CaBeforeEvent{1},1),1); ones(size(Output.CaBeforeEvent{2},1),1) *2; ones(size(Output.CaBeforeEvent{3},1),1) *3],'off');

%check Baseline
nexttile(14)
b = bar([nanmean(Output.CaBaseline{1}),...
    nanmean(Output.CaBaseline{2}),...
    nanmean(Output.CaBaseline{3}),...
    nanmean(Output.CaBaseline{4}),]);
b.FaceColor = 'flat';
b.CData([1],:) = PlotColors{1,1};
b.CData([2],:) = PlotColors{2,1};
b.CData([3],:) = PlotColors{3,1};
b.CData([4],:) = PlotColors{4,1};
ylabel('Mean activity');
title(strcat('Baseline', num2str(Bsl(1,1) -5), 'to', num2str(Bsl(2,1)-5)));
Labels = {'SO','Spindle','SO+spindle'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    Data = Output.CaBaseline{barIndex};
    sorted_data_y = sort(Data, 'ascend');
    sorted_data_x = [];
    [N,edges] = histcounts(Data(:,1),100);
    for iBin = 1: size(N,2)
        if N(1,iBin) >0
            XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
            if N(1,iBin) == 1
                XtmpLocation = 0;
            end
            sorted_data_x = [sorted_data_x XtmpLocation];
        end
    end
    sorted_data_x = sorted_data_x/10+barIndex;
    StandardError(barIndex) = nanstd(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  nanmean(Data);
end
e = errorbar([1:size(StandardError,2)],tmp,StandardError,'.');
clear tmp
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
plotSignificanceBetweenNonParametric(Output.CaBaseline)

ax = nexttile(15,[1 3]);
[p,tbl,stats] = anovan([Output.CaBaseline{1};Output.CaBaseline{2};Output.CaBaseline{3}],...
    {[ones(size(Output.CaBaseline{1},1),1); ones(size(Output.CaBaseline{2},1),1) *2; ones(size(Output.CaBaseline{3},1),1) *3]},'display','off');
VariableNames = tbl(1,:);
ColumnNames = tbl(2:end,1);
tbl = cell2table(tbl);
tbl.Properties.VariableNames = VariableNames;
tbl(1,:) = [];
table2array(tbl);
axisTable(ax, cell2mat(table2array(tbl(1,[3 6 7]))),{'d.f.','F','Prob>F'},ColumnNames(1,1))
axis off
title('One way ANOVA')
p = kruskalwallis([Output.CaDuringEventDuration{1};Output.CaDuringEventDuration{2};Output.CaDuringEventDuration{3}],...
    [ones(size(Output.CaDuringEventDuration{1},1),1); ones(size(Output.CaDuringEventDuration{2},1),1) *2; ones(size(Output.CaDuringEventDuration{3},1),1) *3],'off');
%%
hFig = figure(3);
tiledlayout(4,1)
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 20 500 500])
set(gcf,'color','white')
set(gcf,'name',' %active frames','NumberTitle','off')
for iCluster = 1:3
    nexttile
    [h, edges] = histcounts(Output.ActiveCellsFraction{1,iCluster},0:0.25:10);
    h = h/sum(h);
    bar(edges(2:end),h,'FaceColor',PlotColors{iCluster,1},'EdgeColor',[1 1 1],'FaceAlpha',1)
    ylim([0 1])
    set(gca, 'YScale', 'log')
    ylabel('Probability (%)')
    xlabel('% active cells')
    box(gca,'off')
end
%%
hFig = figure(4);
tiledlayout(4,1)
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 20 500 500])
set(gcf,'color','white')
set(gcf,'name',' %active frames','NumberTitle','off')
for iCluster = 1:3
    nexttile
    [h, edges] = histcounts(Output.EventDurationAll{1,iCluster},0:0.2:2.5);
    h = h/sum(h);
    bar(edges(2:end),h,'FaceColor',PlotColors{iCluster,1},'EdgeColor',[1 1 1],'FaceAlpha',1)
    ylabel('Probability (%)')
    xlabel('EventDuration (s)')
    box(gca,'off')
end