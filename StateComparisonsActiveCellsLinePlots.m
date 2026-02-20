% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
%% Triplet Analysis
close all
clear
UseIntervalBeforeEvent =0;
CellClusterName{1,1} = 'SO-active';
CellClusterName{2,1} = 'Spindle-active';
CellClusterName{3,1} = 'SO+Spindle-active';
MicroaeousalMerge =1;
DirData = '/path/to/data/';
MinEpisodeDuration = 40;
Files = dir(DirData);

for iFile = 3: size(Files,1)
    FileName{iFile-2,1} = Files(iFile,1).name;
end

for iCluster = 1:3
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
    end
    
    CaDataThirds.SWS.Amp = [];
    CaDataThirds.SWS.Freq = [];
    CaDataThirds.REM.Amp = [];
    CaDataThirds.REM.Freq = [];
    CaDataThirds.Wake.Amp = [];
    CaDataThirds.Wake.Freq = [];
    CaDataThirds.SWS.AnimalName = [];
    CaDataThirds.REM.AnimalName = [];
    CaDataThirds.Wake.AnimalName = [];
    CaDataState.SWS.Amp = [];
    CaDataState.SWS.Freq = [];
    CaDataState.REM.Amp = [];
    CaDataState.REM.Freq = [];
    CaDataState.Wake.Amp = [];
    CaDataState.Wake.Freq = [];
    CaDataState.SWS.AnimalName = [];
    CaDataState.REM.AnimalName = [];
    CaDataState.Wake.AnimalName = [];
    CaDataPerCellState.SWS.Amp = [];
    CaDataPerCellState.SWS.Freq = [];
    CaDataPerCellState.REM.Amp = [];
    CaDataPerCellState.REM.Freq = [];
    CaDataPerCellState.Wake.Amp = [];
    CaDataPerCellState.Wake.Freq = [];
    
    for iFile = 1:length(FileName)
        load(strcat(DirData,FileName{iFile,1}),'Events','CaDataAll','SlstAll','FrameRate');
        cfg = [];
        cfg.scoring = SlstAll;
        cfg.scoring_epoch_length = 1;
        cfg.code_NREM = -2;
        cfg.code_REM = -3;
        cfg.code_WAKE =0;
        [Episodes] = SleepEpisodeDetector(cfg);
        
        cfg = [];
        cfg.MinEpisodeDuration = MinEpisodeDuration;
        cfg.MicroarousalMerge = true;
        [Triplets] = ExtractSleepTriplets(Episodes, cfg);
        
        WAKEpisodes = Triplets.WAKEBef;
        WAKThirds = Triplets.WAKEBefThirds;
        NREMEpisodes = Triplets.NREM_before;
        NREMThirds = Triplets.NREM_beforeThirds;
        REMEpisodes = Triplets.REM;
        REMThirds = Triplets.REMThirds;
        NREM_after = Triplets.NREM_after;
        NREM_afterThirds = Triplets.NREM_afterThirds;
        
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
                CurrentEvents = Events.Slo(2:3,SOSpi);
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
        
        if UseIntervalBeforeEvent ==1
            CurrentEvents = [CurrentEvents(1,:)-1;CurrentEvents(1,:)];
        end
        
        [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
        
        %%
        for iState = 1:3
            switch iState
                case 1
                    TMPEpisodes = NREMThirds;
                case 2
                    TMPEpisodes = REMThirds;
                case 3
                    TMPEpisodes = WAKThirds;
            end
            
            TMPCaDataThirdAmp = nan(3,size(setdiff(1:size(TMPEpisodes,2),NoEventNREMEpisodes),2));
            TMPCaDataThirdFre = nan(3,size(setdiff(1:size(TMPEpisodes,2),NoEventNREMEpisodes),2));
            
            for iEpisode = 1:size(TMPEpisodes,2)
                if ~isempty(ActiveCells{1,iEpisode})
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./((TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode))));
                    end
                else
                    TMPCaDataThirdAmp(:,iEpisode) = nan;
                    TMPCaDataThirdFre(:,iEpisode) = nan;
                end
            end
            
            TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
            TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
            TMPCaDataThirdAmp(isnan(TMPCaDataThirdAmp)) = 0;
            
            switch iState
                case 1
                    CaDataThirds.SWS.Amp = [CaDataThirds.SWS.Amp; TMPCaDataThirdAmp'];
                    CaDataThirds.SWS.Freq = [CaDataThirds.SWS.Freq; TMPCaDataThirdFre'];
                    CaDataThirds.SWS.AnimalName = [CaDataThirds.SWS.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdFre,2),1)];
                case 2
                    CaDataThirds.REM.Amp = [CaDataThirds.REM.Amp; TMPCaDataThirdAmp'];
                    CaDataThirds.REM.Freq = [CaDataThirds.REM.Freq; TMPCaDataThirdFre'];
                    CaDataThirds.REM.AnimalName = [CaDataThirds.REM.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdFre,2),1)];
                case 3
                    CaDataThirds.Wake.Amp = [CaDataThirds.Wake.Amp; TMPCaDataThirdAmp'];
                    CaDataThirds.Wake.Freq = [CaDataThirds.Wake.Freq; TMPCaDataThirdFre'];
                    CaDataThirds.Wake.AnimalName = [CaDataThirds.Wake.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdFre,2),1)];
            end
        end
        
        for iState = 1:3
            switch iState
                case 1
                    TMPEpisodes = NREMThirds;
                case 2
                    TMPEpisodes = REMThirds;
                case 3
                    TMPEpisodes = WAKThirds;
            end
            
            TMPCaDataAmp = nan(1,size(setdiff(1:size(TMPEpisodes,2),NoEventNREMEpisodes),2));
            TMPCaDataFre = nan(1,size(setdiff(1:size(TMPEpisodes,2),NoEventNREMEpisodes),2));
            TMPCaDataPerCellAmp = [];
            TMPCaDataPerCellFre = [];
            
            for iEpisode = 1:size(TMPEpisodes,2)
                if ~isempty(ActiveCells{1,iEpisode})
                    TMPCaDataAmp(1,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate)),2));
                    TMPCaDataFre(1,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate))),2)./(TMPEpisodes(4,iEpisode)-TMPEpisodes(1,iEpisode)));
                    TMPCaDataPerCellAmp = [TMPCaDataPerCellAmp; (nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate)),2))];
                    TMPCaDataPerCellFre = [TMPCaDataPerCellFre; (sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate))),2)./(TMPEpisodes(4,iEpisode)-TMPEpisodes(1,iEpisode)))];
                else
                    TMPCaDataAmp(:,iEpisode) = nan;
                    TMPCaDataFre(:,iEpisode) = nan;
                end
            end
            
            TMPCaDataAmp(:,sum(isnan(TMPCaDataFre),1)==size(TMPCaDataFre,1)) = [];
            TMPCaDataFre(:,sum(isnan(TMPCaDataFre),1)==size(TMPCaDataFre,1)) = [];
            TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
            
            switch iState
                case 1
                    CaDataState.SWS.Amp = [CaDataState.SWS.Amp; TMPCaDataAmp'];
                    CaDataState.SWS.Freq = [CaDataState.SWS.Freq; TMPCaDataFre'];
                    CaDataState.SWS.AnimalName = [CaDataState.SWS.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.SWS.Amp = [CaDataPerCellState.SWS.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.SWS.Freq = [CaDataPerCellState.SWS.Freq; TMPCaDataPerCellFre];
                case 2
                    CaDataState.REM.Amp = [CaDataState.REM.Amp; TMPCaDataAmp'];
                    CaDataState.REM.Freq = [CaDataState.REM.Freq; TMPCaDataFre'];
                    CaDataState.REM.AnimalName = [CaDataState.REM.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.REM.Amp = [CaDataPerCellState.REM.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.REM.Freq = [CaDataPerCellState.REM.Freq; TMPCaDataPerCellFre];
                case 3
                    CaDataState.Wake.Amp = [CaDataState.Wake.Amp; TMPCaDataAmp'];
                    CaDataState.Wake.Freq = [CaDataState.Wake.Freq; TMPCaDataFre'];
                    CaDataState.Wake.AnimalName = [CaDataState.Wake.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.Wake.Amp = [CaDataPerCellState.Wake.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.Wake.Freq = [CaDataPerCellState.Wake.Freq; TMPCaDataPerCellFre];
            end
        end
    end
    
    size(CaDataState.REM.Freq)
    size(CaDataState.Wake.Freq)
    
    %% state histograms
    hFig = figure(1);
    if iCluster == 1
        tiledlayout(3,1)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' Across Cells','NumberTitle','off')
    end
    nexttile
    BrainStates{1,1} = 'Wake';
    BrainStates{2,1} = 'SWS';
    BrainStates{3,1} = 'REM';
    Labels = {'Wake', 'SWS', 'REM'};
    PlotColors{1,1} = [0.8500 0.3250 0.0980];
    PlotColors{2,1} = [0 0.4470 0.7410];
    PlotColors{3,1} = [0.3 0.7 0.1];
    
    for iState = 1:3
        Data = CaDataPerCellState.(BrainStates{iState,1}).Freq;
        [h, edges] = histcounts(log(Data(Data>0)),[-7:0.5:0]);
        h = h/sum(h);
        bar(edges(2:end),h,'FaceColor',PlotColors{iState,1},'EdgeColor',[1 1 1],'FaceAlpha',0.2)
        ft = fittype( 'gauss1' );
        fitresult = fit(edges(2:end)', h', ft);
        hold on
        plt = plot(fitresult);
        plt.Color= PlotColors{iState,1};
        plt.LineWidth = 1.5;
    end
    
    xlim([-7 0])
    ylabel('Probability (%)')
    xlabel('log Frequency of Ca^{2+} transients (Hz)')
    legend('Wake','','SWS','','REM')
    legend boxoff
    box(gca,'off')
    
    hFig = figure(2);
    if iCluster == 1
        tiledlayout(3,1)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' Across Episodes','NumberTitle','off')
    end
    
    nexttile
    BrainStates{1,1} = 'Wake';
    BrainStates{2,1} = 'SWS';
    BrainStates{3,1} = 'REM';
    Labels = {'Wake', 'SWS', 'REM'};
    PlotColors{1,1} = [0.3 0.7 0.1];
    PlotColors{2,1} = [0.8500 0.3250 0.0980];
    PlotColors{3,1} = [0 0.4470 0.7410];
    
    for iState = 1:3
        Data = CaDataState.(BrainStates{iState,1}).Freq;
        [h, edges] = histcounts(Data,0:0.005:0.08);
        h = h/sum(h);
        bar(edges(2:end),h,'FaceColor',PlotColors{iState,1},'EdgeColor',[1 1 1],'FaceAlpha',0.2)
        ft = fittype( 'smoothingspline' );
        fitresult = fit(edges(2:end)', h', ft);
        hold on
        plt = plot(fitresult);
        plt.Color= PlotColors{iState,1};
        plt.LineWidth = 1.5;
    end
    
    ylim([0 0.3])
    xlim([0 0.1])
    ylabel('Probability (%)')
    xlabel('log Frequency of Ca^{2+} transients (Hz)')
    legend(' ','Wake','','SWS','','REM')
    legend boxoff
    box(gca,'off')
    
    hFig = figure(33);
    if iCluster == 1
        tiledlayout(6,1)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 200 1000])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
    end
    
    BrainStates{1,1} = 'Wake';
    BrainStates{2,1} = 'SWS';
    BrainStates{3,1} = 'REM';
    Labels = {'Wake', 'SWS', 'REM'};
    PlotColors{1,1} = [0.8500 0.3250 0.0980];
    PlotColors{2,1} = [0 0.4470 0.7410];
    PlotColors{3,1} = [0.3 0.7 0.1];
    
    nexttile
    Data = [];
    for iState = 1:3
        Data(:,iState) = CaDataState.(BrainStates{iState,1}).Freq;
    end
    Data(sum(isnan(Data),2)==3,:) = [];
    Data(isnan(Data))=0;
    
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColors{iCluster,1});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
    ylim([0 0.04])
    xlim([0.5 3.5])
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    plotSignificanceWithinNonParametric(Data)
    
    nexttile
    Data = [];
    for iState = 1:3
        Data(:,iState) = CaDataState.(BrainStates{iState,1}).Amp*100;
    end
    Data(sum(isnan(Data),2)==3,:) = [];
    Data(isnan(Data))=0;
    
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColors{iCluster,1});
    title('Active Cells');
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    ylim([0 150])
    xlim([0.5 3.5])
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    plotSignificanceWithinNonParametric(Data)
    
    %% state thirds
    hFig = figure(44);
    if iCluster == 1
        tiledlayout(6,3)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 1000])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
    end
    BrainStates{1,1} = 'Wake';
    BrainStates{2,1} = 'SWS';
    BrainStates{3,1} = 'REM';
    Labels = {'1^{st}', '2^{nd}', '3^{rd}'};
    PlotColors{1,1} = [0.8500 0.3250 0.0980];
    PlotColors{2,1} = [0 0.4470 0.7410];
    PlotColors{3,1} = [0.3 0.7 0.1];
    
    for iState = 1:3
        nexttile
        Data = CaDataThirds.(BrainStates{iState,1}).Freq;
        Data(sum(isnan(Data),2)==3,:) = [];
        Data(isnan(Data))=0;
        
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColors{iCluster,1});
        title(BrainStates{iState,1});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 0.04])
        xlim([0.5 3.5])
        xlabel('Third')
        box(gca,'off')
        plotSignificanceWithinNonParametric(Data)
    end
    
    for iState = 1:3
        Data = CaDataThirds.(BrainStates{iState,1}).Amp*100;
        nexttile       
        Data(sum(isnan(Data),2)==3,:) = [];
        Data(isnan(Data))=0;
        
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColors{iCluster,1});
        title(BrainStates{iState,1});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
        ylim([0 150])
        xlim([0.5 3.5])
        xlabel('Third')
        plotSignificanceWithinNonParametric(Data)
    end
    box(gca,'off')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    
    BrainStates{1,1} = 'SWS';
    BrainStates{2,1} = 'REM';
    BrainStates{3,1} = 'Wake';
    
    %prepare data within states
    StateCond{iCluster,1} = [];
    Third{iCluster,1} = [];
    DataForLMMFreq{iCluster,1} = [];
    DataForLMMAmp{iCluster,1} = [];
    AnimalName{iCluster,1} = [];
    Cluster{iCluster,1} = [];
    
    for iState = 1:3
        DataForLMMFreq{iCluster,1} = [DataForLMMFreq{iCluster,1}; [CaDataThirds.(BrainStates{iState,1}).Freq(:,1);...
            CaDataThirds.(BrainStates{iState,1}).Freq(:,2);...
            CaDataThirds.(BrainStates{iState,1}).Freq(:,3)]];
        DataForLMMAmp{iCluster,1} = [DataForLMMAmp{iCluster,1}; CaDataThirds.(BrainStates{iState,1}).Amp(:,1);...
            CaDataThirds.(BrainStates{iState,1}).Amp(:,2);...
            CaDataThirds.(BrainStates{iState,1}).Amp(:,3)];
        StateCond{iCluster,1} = [StateCond{iCluster,1}; repmat(iState,size(CaDataThirds.(BrainStates{iState,1}).Freq,1)*3,1)];
        Third{iCluster,1} = [Third{iCluster,1}; ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1);...
            ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1) *2;...
            ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1)*3];
        AnimalName{iCluster,1} = [AnimalName{iCluster,1}; repmat(CaDataThirds.(BrainStates{iState,1}).AnimalName,3,1)];
    end
    Cluster{iCluster,1} = [repmat(iCluster,length(AnimalName{iCluster,1}),1)];
    
    %prepare data across states
    StateCond{iCluster,2} = [];
    Third{iCluster,2} = [];
    DataForLMMFreq{iCluster,2} = [];
    DataForLMMAmp{iCluster,2} = [];
    AnimalName{iCluster,2} = [];
    Cluster{iCluster,2} = [];
    
    for iState = 1:3
        DataForLMMFreq{iCluster,2} = [DataForLMMFreq{iCluster,2}; [CaDataThirds.(BrainStates{iState,1}).Freq(:,1);...
            CaDataThirds.(BrainStates{iState,1}).Freq(:,2);...
            CaDataThirds.(BrainStates{iState,1}).Freq(:,3)]];
        DataForLMMAmp{iCluster,2} = [DataForLMMAmp{iCluster,2}; CaDataThirds.(BrainStates{iState,1}).Amp(:,1);...
            CaDataThirds.(BrainStates{iState,1}).Amp(:,2);...
            CaDataThirds.(BrainStates{iState,1}).Amp(:,3)];
        StateCond{iCluster,2} = [StateCond{iCluster,2}; repmat(iState,size(CaDataThirds.(BrainStates{iState,1}).Freq,1)*3,1)];
        Third{iCluster,2} = [Third{iCluster,2}; ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1);...
            ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1) *2;...
            ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1)*3];
        AnimalName{iCluster,2} = [AnimalName{iCluster,2}; repmat(CaDataThirds.(BrainStates{iState,1}).AnimalName,3,1)];
    end
    Cluster{iCluster,2} = [repmat(iCluster,length(AnimalName{iCluster,2}),1)];
    
    StateCond{iCluster,3} = [];
    Third{iCluster,3} = [];
    DataForLMMFreq{iCluster,3} = [];
    DataForLMMAmp{iCluster,3} = [];
    AnimalName{iCluster,3} = [];
    Cluster{iCluster,3} = [];
    
    for iState = 1:3
        DataForLMMFreq{iCluster,3} = [DataForLMMFreq{iCluster,3}; [CaDataState.(BrainStates{iState,1}).Freq(:,1)]];
        DataForLMMAmp{iCluster,3} = [DataForLMMAmp{iCluster,3}; CaDataState.(BrainStates{iState,1}).Amp(:,1)];
        StateCond{iCluster,3} = [StateCond{iCluster,3}; repmat(iState,size(CaDataState.(BrainStates{iState,1}).Freq,1),1)];
        Third{iCluster,3} = [Third{iCluster,3}; ones(size(CaDataState.(BrainStates{iState,1}).Freq,1),1)];
        AnimalName{iCluster,3} = [AnimalName{iCluster,3}; repmat(CaDataState.(BrainStates{iState,1}).AnimalName,1,1)];
    end
    Cluster{iCluster,3} = [repmat(iCluster,length(AnimalName{iCluster,3}),1)];
end

%across states
T = [];
for iCluster = 1:3
    T = [T;array2table([AnimalName{iCluster,3}, DataForLMMFreq{iCluster,3}, StateCond{iCluster,3}, Third{iCluster,3},Cluster{iCluster,3} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds','Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);
T.Cluster = categorical(T.Cluster);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Cluster + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test interaction
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect State
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect Cluster
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

T = [];
for iCluster = 1:3
    T = [T;array2table([AnimalName{iCluster,3}, DataForLMMAmp{iCluster,3}, StateCond{iCluster,3}, Third{iCluster,3},Cluster{iCluster,3} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds','Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);
T.Cluster = categorical(T.Cluster);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Cluster + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test interaction
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect State
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect Cluster
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

%across states with triplets
T = [];
for iCluster = 1:3
    T = [T;array2table([AnimalName{iCluster,2}, DataForLMMFreq{iCluster,2}, StateCond{iCluster,2}, Third{iCluster,2},Cluster{iCluster,2} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds','Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);
T.Cluster = categorical(T.Cluster);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Cluster + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test interaction
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect State
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); %test Main Effect Cluster
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

%within states
T = [];
for iCluster = 1:3
    T = [T;array2table([AnimalName{iCluster,1}, DataForLMMFreq{iCluster,1}, StateCond{iCluster,1}, Third{iCluster,1},Cluster{iCluster,1} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds','Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);
T.Cluster = categorical(T.Cluster);

for iState = 1:3
    T_subset = T(double(T.State) == iState, :);
    lm1Model = fitlme(T_subset,strcat('Freq',' ~ Thirds *Cluster + (1|Animal)'));
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Thirds +Cluster + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test interaction
    if results.pValue(end) <= 0.05
        PValues(1,1) = 1 - results.pValue(end);
    end
    
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Cluster + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test Main Effect Thirds
    if results.pValue(end) <= 0.05
        PValuesDiet(1,1) = 1 - results.pValue(end);
    end
    
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Thirds + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test Main Effect Cluster
    if results.pValue(end) <= 0.05
        PValuesGlucose(1,1) = 1 - results.pValue(end);
    end
end

%within states
T = [];
for iCluster = 1:3
    T = [T;array2table([AnimalName{iCluster,1}, DataForLMMAmp{iCluster,1}, StateCond{iCluster,1}, Third{iCluster,1},Cluster{iCluster,1} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds','Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);
T.Cluster = categorical(T.Cluster);

for iState = 1:3
    T_subset = T(double(T.State) == iState, :);
    lm1Model = fitlme(T_subset,strcat('Freq',' ~ Thirds *Cluster + (1|Animal)'));
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Thirds +Cluster + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test interaction
    if results.pValue(end) <= 0.05
        PValues(1,1) = 1 - results.pValue(end);
    end
    
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Cluster + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test Main Effect Thirds
    if results.pValue(end) <= 0.05
        PValuesDiet(1,1) = 1 - results.pValue(end);
    end
    
    lm2Model = fitlme(T_subset,strcat('Freq',' ~ Thirds + (1|Animal)'));
    results = compare(lm2Model, lm1Model); %test Main Effect Cluster
    if results.pValue(end) <= 0.05
        PValuesGlucose(1,1) = 1 - results.pValue(end);
    end
end

% Full model with all interactions
lmFull = fitlme(T, 'Freq ~ State * Thirds * Cluster + (1|Animal)');

% 1) Test 3-way interaction
lmReduced = fitlme(T, 'Freq ~ State * Thirds + State * Cluster + Thirds * Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.Interaction_3way = results.pValue(end);

% If not significant, drop 3-way interaction
if results.pValue(end) > 0.05
    lmFull = lmReduced;
end

% 2) Test State × Thirds interaction
lmReduced = fitlme(T, 'Freq ~ State + Thirds + Cluster + State * Cluster + Thirds * Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.State_Thirds = results.pValue(end);
if results.pValue(end) > 0.05
    lmFull = lmReduced;
end

% 3) Test State × Cluster interaction
lmReduced = fitlme(T, 'Freq ~ State + Thirds + Cluster + State:Thirds + Thirds:Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.State_Cluster = results.pValue(end);
if results.pValue(end) > 0.05
    lmFull = lmReduced;
end

% 4) Test Thirds × Cluster interaction
lmReduced = fitlme(T, 'Freq ~ State + Thirds + Cluster + State * Thirds + State * Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.Thirds_Cluster = results.pValue(end);
if results.pValue(end) > 0.05
    lmFull = lmReduced;
end

% 5) Test main effect of State
lmReduced = fitlme(T, 'Freq ~ Thirds + Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.State_Main = results.pValue(end);

% 6) Test main effect of Thirds
lmReduced = fitlme(T, 'Freq ~ State + Cluster + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.Thirds_Main = results.pValue(end);

% 7) Test main effect of Cluster
lmReduced = fitlme(T, 'Freq ~ State + Thirds + (1|Animal)');
results = compare(lmReduced, lmFull);
PValues.Cluster_Main = results.pValue(end);