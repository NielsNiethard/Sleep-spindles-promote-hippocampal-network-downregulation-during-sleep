% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
%% Triplet Analysis
close all
clear
UseIntervalBeforeEvent =1;
IncludePreceedingWake =0;
MicroaeousalMerge =1;
DirData = '/path/to/data/';
MinEpisodeDuration = 40;
Files = dir(DirData);

for iFile = 3: size(Files,1)
    FileName{iFile-2,1} = Files(iFile,1).name;
end

PlotColor{1,1} = [1 0 0];
PlotColor{2,1} = [0 1 0];
PlotColor{3,1} = [0 0 1];
CellClusterName{1,1} = 'SOs';
CellClusterName{2,1} = 'Spindles';
CellClusterName{3,1} = 'SO+Spindles';

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
    
    EpisodeNREMNDuration.REM = [];
    EpisodeNREMN1Duration.REM = [];
    EpisodeNREMNDuration.WAK = [];
    EpisodeNREMN1Duration.WAK = [];
    EventsDensity.WAK = [];
    EventsPerEpisode.REM = [];
    EventsPerEpisode.WAK = [];
    EpisodeDistance.REM = [];
    EpisodeDistance.WAK = [];
    CaDataThirdAmp = [];
    CaDataThirdFre =[];
    CaDataThirds.SWSREMSWS.Amp = [];
    CaDataThirds.SWSREMSWS.Freq = [];
    CaData.SWSREMSWS.Amp = [];
    CaData.SWSREMSWS.Freq = [];
    CaDataThirds.SWSWAKSWS.Amp = [];
    CaDataThirds.SWSWAKSWS.Freq = [];
    CaData.SWSWAKSWS.Amp = [];
    CaData.SWSWAKSWS.Freq = [];
    
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
        
        NREMEpisodes = Triplets.NREM_before;
        NREMThirds = Triplets.NREM_beforeThirds;
        REMEpisodes = Triplets.REM;
        REMThirds = Triplets.REMThirds;
        NREMEpisodesAfter = Triplets.NREM_after;
        NREMEpisodesThirdsAfter = Triplets.NREM_afterThirds;
        
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
        
        [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
        
        %%
        cfg = [];
        cfg.MinEpisodeDuration = MinEpisodeDuration;
        cfg.MicroarousalMerge = true;
        cfg.Mode = 'WAKE';
        [Triplets] = ExtractSleepTriplets(Episodes, cfg);
        
        WAKEpisodes = Triplets.WAKE;
        WAKThirds = Triplets.WAKEThirds;
        NREMEpisodes = Triplets.NREM_before;
        NREMThirds = Triplets.NREM_beforeThirds;
        REMEpisodes = Triplets.REM;
        REMThirds = Triplets.REMThirds;
        NREMEpisodesAfter = Triplets.NREM_after;
        NREMEpisodesThirdsAfter = Triplets.NREM_afterThirds;
        
        TMPCaDataThirdAmp = nan(9,size(NREMThirds,2));
        TMPCaDataThirdFre = nan(9,size(NREMThirds,2));
        TMPCaDataAmp = nan(1,size(NREMThirds,2));
        TMPCaDataFre = nan(1,size(NREMThirds,2));
        WAKEpisodeDistance = nan(1,size(NREMThirds,2));
        WAKEpisodeNREMNDuration = nan(1,size(NREMThirds,2));
        WAKEpisodeNREMN1Duration = nan(1,size(NREMThirds,2));
        
        for iTriplet = 1
            for iEpisode = 1:size(NREMThirds,2)
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMThirds;
                    case 2
                        TMPEpisodes = WAKThirds;
                    case 3
                        TMPEpisodes = NREMEpisodesThirdsAfter;
                end
                
                if ~isnan(nanmean(nanmean(nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2))))
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(:,round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                    end
                    
                    switch iTriplet
                        case 1
                            TMPEpisodes = NREMEpisodes;
                        case 2
                            TMPEpisodes = WAKEpisodes;
                        case 3
                            TMPEpisodes = NREMEpisodesAfter;
                    end
                    
                    TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2));
                    TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                    
                    if iTriplet ==1
                        WAKEpisodeDistance(1,iEpisode) = TMPEpisodes(2,iEpisode) - TMPEpisodes(1,iEpisode);
                        WAKEpisodeNREMNDuration(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    end
                end
            end
        end
        
        WAKEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        NumberOfEventsPerEpisode(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataFre(:,sum(isnan(TMPCaDataFre),1)==size(TMPCaDataFre,1)) = [];
        
        CaDataThirds.SWSWAKSWS.Amp = [CaDataThirds.SWSWAKSWS.Amp; TMPCaDataThirdAmp'];
        CaDataThirds.SWSWAKSWS.Freq = [CaDataThirds.SWSWAKSWS.Freq; TMPCaDataThirdFre'];
        EpisodeNREMNDuration.WAK = [EpisodeNREMNDuration.WAK , WAKEpisodeNREMNDuration(:,:)];
        EpisodeNREMN1Duration.WAK = [EpisodeNREMN1Duration.WAK, WAKEpisodeNREMN1Duration(:,:)];
        CaData.SWSWAKSWS.Amp = [CaData.SWSWAKSWS.Amp; TMPCaDataAmp(:,:)'];
        CaData.SWSWAKSWS.Freq = [CaData.SWSWAKSWS.Freq; TMPCaDataFre(:,:)'];
        EpisodeDistance.WAK = [EpisodeDistance.WAK, WAKEpisodeDistance(1,:)];
        EventsPerEpisode.WAK = [EventsPerEpisode.WAK; NumberOfEventsPerEpisode(1,:)'];
        EventsDensity.WAK = [EventsDensity.WAK; NumberOfEventsPerEpisode(1,:)'./ ...
            WAKEpisodeNREMNDuration'];
    end
    
    %%
    %Mean across episodes
    hFig = figure(200);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [100 20 1600 1000])
        set(gcf,'color','white')
        set(gcf,'name',' Summary All Cells','NumberTitle','off')
    end
    nexttile
    
    %correlations
    CorrType = 'Spearman';
    Color = 'k';
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensity.WAK,CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Density ',CellClusterName{iCluster,1}));
    ylabel ('Frequency')
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensity.WAK,CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Density ',CellClusterName{iCluster,1}));
    ylabel ('Amplitude')
    title ('Amp')
    DataDensity{:,iCluster} = EventsDensity.WAK;
end

hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [100 100 200 200])
set(gcf,'color','white')
set(gcf,'name',' Density','NumberTitle','off')
b = bar([nanmean(DataDensity{1,1}),...
    nanmean(DataDensity{1,2}),...
    nanmean(DataDensity{1,3}),]);
b.FaceColor = 'flat';
b.CData([1],:) = PlotColor{1,1};
b.CData([2],:) = PlotColor{2,1};
b.CData([3],:) = PlotColor{3,1};
ylabel('Density (Hz)');
title('EventDensity');
Labels = {'SO','Spindle','SO+spindle'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    Data = DataDensity{1,barIndex};
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
plotSignificanceBetweenTTest(DataDensity)