% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
%% Triplet Analysis
close all
clear
UseIntervalBeforeEvent =0;
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
PlotColor{4,1} = [0.3 0.3 0.3];
CellClusterName{1,1} = 'SOs';
CellClusterName{2,1} = 'Spindles';
CellClusterName{3,1} = 'SO+Spindles';
CellClusterName{4,1} = 'Control';

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
    end
    
    EpisodeNREMNDuration.REM = [];
    EpisodeNREMN1Duration.REM = [];
    EpisodeNREMNDuration.WAK = [];
    EpisodeNREMN1Duration.WAK = [];
    EventsPerEpisode.REM = [];
    EventsPerEpisode.WAK = [];
    EventsDensityPerEpisode.REM = [];
    EventsDensityPerEpisode.WAK = [];
    EpisodeDistance.REM = [];
    EpisodeDistance.WAK = [];
    CaDataThirdAmp = [];
    CaDataThirdFre =[];
    CaDataThirds.SWSREMSWS.Amp = [];
    CaDataThirds.SWSREMSWS.Freq = [];
    CaData.SWSREMSWS.Amp = [];
    CaData.SWSREMSWS.Freq = [];
    CaData.SWSREMSWS.AnimalNameFreq = [];
    CaData.SWSREMSWS.AnimalNameAmp = [];
    CaDataThirds.SWSWAKSWS.Amp = [];
    CaDataThirds.SWSWAKSWS.Freq = [];
    CaDataThirds.SWSREMSWS.AnimalNameFreq = [];
    CaDataThirds.SWSREMSWS.AnimalNameAmp = [];
    CaDataThirds.SWSWAKSWS.AnimalNameFreq = [];
    CaDataThirds.SWSWAKSWS.AnimalNameAmp = [];
    CaData.SWSWAKSWS.Amp = [];
    CaData.SWSWAKSWS.Freq = [];
    CaData.SWSWAKSWS.AnimalNameFreq = [];
    CaData.SWSWAKSWS.AnimalNameAmp = [];
    
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
        
        if iCluster <= 3
            [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
        else
            CurrentEvents = [Events.Spi, [Events.Slo(4,:);Events.Slo(3,:)]];
            [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
            ActiveCells=InactiveCells;
        end
        
        %%
        TMPCaDataThirdAmp = nan(9,size(REMThirds,2));
        TMPCaDataThirdFre = nan(9,size(REMThirds,2));
        TMPCaDataAmp = nan(3,size(REMThirds,2));
        TMPCaDataFre = nan(3,size(REMThirds,2));
        REMEpisodeDistance = nan(1,size(REMThirds,2));
        REMEpisodeNREMNDuration = nan(1,size(REMThirds,2));
        REMEpisodeNREMN1Duration = nan(1,size(REMThirds,2));
        REMEventDensity = nan(1,size(REMThirds,2));
        
        for iTriplet = 1:3
            for iEpisode = 1:size(REMThirds,2)
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMThirds;
                    case 2
                        TMPEpisodes = REMThirds;
                    case 3
                        TMPEpisodes = NREMEpisodesThirdsAfter;
                end
                
                for iThird = 1:3
                    TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                    TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./((TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode))));
                end
                
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMEpisodes;
                    case 2
                        TMPEpisodes = REMEpisodes;
                    case 3
                        TMPEpisodes = NREMEpisodesAfter;
                end
                
                TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2));
                TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                
                if iTriplet ==2
                    REMEpisodeDistance(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    REMEpisodeNREMNDuration(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    REMEpisodeNREMN1Duration(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodesAfter(1,iEpisode);
                    REMEventDensity(1,iEpisode) = EventDensity(1,iEpisode);
                end
            end
        end
        
        REMEventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        NumberOfEventsPerEpisode(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdFre(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataThirdAmp(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataFre(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
        
        EpisodeNREMNDuration.REM = [EpisodeNREMNDuration.REM , REMEpisodeNREMNDuration];
        EpisodeNREMN1Duration.REM = [EpisodeNREMN1Duration.REM, REMEpisodeNREMN1Duration];
        EpisodeDistance.REM = [EpisodeDistance.REM,REMEpisodeDistance];
        EventsPerEpisode.REM = [EventsPerEpisode.REM; NumberOfEventsPerEpisode'];
        EventsDensityPerEpisode.REM = [EventsDensityPerEpisode.REM; REMEventDensity'];
        
        CaDataThirds.SWSREMSWS.Amp = [CaDataThirds.SWSREMSWS.Amp; TMPCaDataThirdAmp'];
        CaDataThirds.SWSREMSWS.Freq = [CaDataThirds.SWSREMSWS.Freq; TMPCaDataThirdFre'];
        CaDataThirds.SWSREMSWS.AnimalNameFreq = [CaDataThirds.SWSREMSWS.AnimalNameFreq; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdFre,2),1)];
        CaDataThirds.SWSREMSWS.AnimalNameAmp = [CaDataThirds.SWSREMSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdAmp,2),1)];
                    
        CaData.SWSREMSWS.Amp = [CaData.SWSREMSWS.Amp; TMPCaDataAmp'];
        CaData.SWSREMSWS.Freq = [CaData.SWSREMSWS.Freq; TMPCaDataFre'];
        CaData.SWSREMSWS.AnimalNameFreq = [CaData.SWSREMSWS.AnimalNameFreq; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
        CaData.SWSREMSWS.AnimalNameAmp = [CaData.SWSREMSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataAmp,2),1)];
        
        %%
        %SWS Wake SWS triplets
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
        
        [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
        
        %%
        TMPCaDataThirdAmp = nan(9,size(WAKThirds,2));
        TMPCaDataThirdFre = nan(9,size(WAKThirds,2));
        TMPCaDataAmp = nan(3,size(WAKThirds,2));
        TMPCaDataFre = nan(3,size(WAKThirds,2));
        WAKEpisodeDistance = nan(1,size(WAKThirds,2));
        WAKEpisodeNREMNDuration = nan(1,size(WAKThirds,2));
        WAKEpisodeNREMN1Duration = nan(1,size(WAKThirds,2));
        WakeEventDensity = nan(1,size(WAKThirds,2));
        
        for iTriplet = 1:3
            for iEpisode = 1:size(WAKThirds,2)
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMThirds;
                    case 2
                        TMPEpisodes = WAKThirds;
                    case 3
                        TMPEpisodes = NREMEpisodesThirdsAfter;
                end
                
                for iThird = 1:3
                    TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                    TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./((TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode))));
                end
                
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMEpisodes;
                    case 2
                        TMPEpisodes = WAKEpisodes;
                    case 3
                        TMPEpisodes = NREMEpisodesAfter;
                end
                
                TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2));
                TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                
                if iTriplet ==2
                    WAKEpisodeDistance(1,iEpisode) = TMPEpisodes(2,iEpisode) - TMPEpisodes(1,iEpisode);
                    WakeEventDensity(1,iEpisode) = EventDensity(1,iEpisode);
                    WAKEpisodeNREMNDuration(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    WAKEpisodeNREMN1Duration(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodesAfter(1,iEpisode);
                end
            end
        end
        
        WakeEventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        NumberOfEventsPerEpisode(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdFre(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataThirdAmp(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataFre(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
        
        %match WAk and REM episode durations
        D=abs(WAKEpisodeDistance(:)-REMEpisodeDistance(:).');
        MatchedWakEpisodes=sortrows(matchpairs(D,max(D(:))),1);
        
        CaDataThirds.SWSWAKSWS.Amp = [CaDataThirds.SWSWAKSWS.Amp; TMPCaDataThirdAmp(:,MatchedWakEpisodes(:,1))'];
        CaDataThirds.SWSWAKSWS.Freq = [CaDataThirds.SWSWAKSWS.Freq; TMPCaDataThirdFre(:,MatchedWakEpisodes(:,1))'];
        CaDataThirds.SWSWAKSWS.AnimalNameAmp = [CaDataThirds.SWSWAKSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdAmp(:,MatchedWakEpisodes(:,1)),2),1)];
        CaDataThirds.SWSWAKSWS.AnimalNameFreq = [CaDataThirds.SWSWAKSWS.AnimalNameFreq; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataThirdFre(:,MatchedWakEpisodes(:,1)),2),1)];
        
        EpisodeNREMNDuration.WAK = [EpisodeNREMNDuration.WAK , WAKEpisodeNREMNDuration(:,MatchedWakEpisodes(:,1))];
        EpisodeNREMN1Duration.WAK = [EpisodeNREMN1Duration.WAK, WAKEpisodeNREMN1Duration(:,MatchedWakEpisodes(:,1))];
        CaData.SWSWAKSWS.Amp = [CaData.SWSWAKSWS.Amp; TMPCaDataAmp(:,MatchedWakEpisodes(:,1))'];
        CaData.SWSWAKSWS.Freq = [CaData.SWSWAKSWS.Freq; TMPCaDataFre(:,MatchedWakEpisodes(:,1))'];
        CaData.SWSWAKSWS.AnimalNameFreq = [CaData.SWSWAKSWS.AnimalNameFreq; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre(:,MatchedWakEpisodes(:,1)),2),1)];
        CaData.SWSWAKSWS.AnimalNameAmp = [CaData.SWSWAKSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataAmp(:,MatchedWakEpisodes(:,1)),2),1)];
        EpisodeDistance.WAK = [EpisodeDistance.WAK, WAKEpisodeDistance(1,MatchedWakEpisodes(:,1))];
        EventsPerEpisode.WAK = [EventsPerEpisode.WAK; NumberOfEventsPerEpisode(1,MatchedWakEpisodes(:,1))'];
        EventsDensityPerEpisode.WAK = [EventsDensityPerEpisode.WAK; WakeEventDensity(1,MatchedWakEpisodes(:,1))'];
    end
    
    %%
    %Mean across episodes
    hFig = figure(200);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [100 20 1600 1000])
        set(gcf,'color','white')
        set(gcf,'name',' Summary All Cells','NumberTitle','off')
        tiledlayout(6,10);
        Labels = {'SWS_{n}', 'REM', 'SWS_{n+1}'};
        PlotColors{1,1} = [0 0.4470 0.7410];
        PlotColors{2,1} = [0.3 0.7 0.1];
        PlotColors{3,1} = [0 0.4470 0.7410];
        Data = CaData.SWSREMSWS.Freq;
        Data(sum(isnan(Data),2)==3,:) =[];
        Data(isnan(Data)) = 0;
        
        nexttile(20*(iCluster-1)+1)
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 0.02])
        box(gca,'off')
        
        [pValue, tbl, stats] = friedman(Data, 1, 'off');
        if pValue < 0.05
            plotSignificanceWithinNonParametric(Data)
        end
        
        %SWS WAKE SWS Triplets
        Labels = {'SWS_{n}', 'WAKE', 'SWS_{n+1}'};
        PlotColors{1,1} = [0 0.4470 0.7410];
        PlotColors{2,1} = [0.8500 0.3250 0.0980];
        PlotColors{3,1} = [0 0.4470 0.7410];
        Data = CaData.SWSWAKSWS.Freq;
        Data(sum(isnan(Data),2)==9,:) =[];
        Data(isnan(Data)) = 0;
        
        nexttile(20*(iCluster-1)+2)
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
        title({'SWS_{n}  WAKE  SWS_{n+1}'});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 0.02])
        box(gca,'off')
        
        %stats
        % Clean up Data1
        Data1 = CaData.SWSREMSWS.Freq;
        Data1(isnan(Data1)) = 0;  
        Data1 = Data1(~any(isnan(Data1), 2), :);  
        % Clean up Data2
        Data2 = CaData.SWSWAKSWS.Freq;
        Data2(isnan(Data2)) = 0;  
        Data2 = Data2(~any(isnan(Data2), 2), :);  
        
        T = array2table([Data1, Data2]);
        T.Properties.VariableNames = {'SWSnREM', 'REM','SWSnnREM', 'SWSnWAK', 'WAK','SWSnnWAK' };
        withinDesign = table([1 1 1 2 2 2]',[1 2 3 1 2 3]','VariableNames',{'State','Triplet'});
        withinDesign.State = categorical(withinDesign.State);
        withinDesign.Triplet = categorical(withinDesign.Triplet);
        
        rm = fitrm(T, 'SWSnREM-SWSnnWAK ~ 1', 'WithinDesign', withinDesign);
        AT = ranova(rm, 'WithinModel', 'State*Triplet');
        disp(anovaTable(AT, 'Value'));
        multcompare(rm, 'State','By','Triplet');
        
        %prepare ANOVA output
        ax= nexttile(20*(iCluster-1)+3,[2,4]);
        axisTable(ax, AT{:,[2 4 5]},{'DF','F','P'},AT.Properties.RowNames')
        axis off
        title('Frequency')
        
        %%
        Data = CaData.SWSREMSWS.Amp*100;
        Data(sum(isnan(Data),2)==9,:) =[];
        Data(isnan(Data)) = 0;
        
        nexttile(20*(iCluster-1)+11)
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
        Labels = {'SWS_{n}', 'REM', 'SWS_{n+1}'};
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 160])
        box(gca,'off')
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
        
        [pValue, tbl, stats] = friedman(Data, 1, 'off');
        if pValue < 0.05
            plotSignificanceWithinNonParametric(Data)
        end
        
        Data = CaData.SWSWAKSWS.Amp*100;
        Data(sum(isnan(Data),2)==9,:) =[];
        Data(isnan(Data)) = 0;
        
        nexttile(20*(iCluster-1)+12)
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
        title({'SWS_{n}  WAKE  SWS_{n+1}'});
        Labels = {'SWS_{n}', 'WAKE', 'SWS_{n+1}'};
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        box(gca,'off')
        
        [pValue, tbl, stats] = friedman(Data, 1, 'off');
        if pValue < 0.05
            plotSignificanceWithinNonParametric(Data)
        end
        ylim([0 160])
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
        
        %stats
        Data1 = CaData.SWSREMSWS.Amp;
        Data1(isnan(Data1)) = 0;  
        Data1 = Data1(~any(isnan(Data1), 2), :);  
        % Clean up Data2
        Data2 = CaData.SWSWAKSWS.Amp;
        Data2(isnan(Data2)) = 0;  
        Data2 = Data2(~any(isnan(Data2), 2), :);  
        
        % Make Data1 and Data2 the same length by adding means to Data1
        if size(Data2, 1) > size(Data1, 1)
            Data1 = [Data1; repmat(mean(Data1, 1), size(Data2, 1) - size(Data1, 1), 1)];
        end
        
        T = array2table([Data1, Data2]);
        T.Properties.VariableNames = {'SWSnREM', 'REM','SWSnnREM', 'SWSnWAK', 'WAK','SWSnnWAK' };
        withinDesign = table([1 1 1 2 2 2]',[1 2 3 1 2 3]','VariableNames',{'State','Triplet'});
        withinDesign.State = categorical(withinDesign.State);
        withinDesign.Triplet = categorical(withinDesign.Triplet);
        
        rm = fitrm(T, 'SWSnREM-SWSnnWAK ~ 1', 'WithinDesign', withinDesign);
        AT = ranova(rm, 'WithinModel', 'State*Triplet');
        disp(anovaTable(AT, 'Value'));
        multcompare(rm, 'State','By','Triplet');
        
        %prepare ANOVA output
        ax= nexttile(20*(iCluster-1)+7,[2,4]);
        axisTable(ax, AT{:,[2 4 5]},{'DF','F','P'},AT.Properties.RowNames')
        axis off
        title('Amplitude')
        
        DiffREMTriplets.Freq{iCluster} = CaData.SWSREMSWS.Freq(:,3)- CaData.SWSREMSWS.Freq(:,1);
        DiffREMTriplets.Amp{iCluster} = CaData.SWSREMSWS.Amp(:,3)- CaData.SWSREMSWS.Amp(:,1);
        DiffWAKTriplets.Freq{iCluster} = CaData.SWSWAKSWS.Freq(:,3)- CaData.SWSWAKSWS.Freq(:,1);
        DiffWAKTriplets.Amp{iCluster} = CaData.SWSWAKSWS.Amp(:,3)- CaData.SWSWAKSWS.Amp(:,1);
        BrainStates{1,1} ='SWSREMSWS';
        BrainStates{2,1} = 'SWSWAKSWS';
        
        %prepare data across Triplets
        StateCond{iCluster,1} = [];
        Third{iCluster,1} = [];
        DataForLMMFreq{iCluster,1} = [];
        DataForLMMAmp{iCluster,1} = [];
        AnimalName{iCluster,1} = [];
        Cluster{iCluster,1} = [];
        
        for iState = 1:2
            DataForLMMFreq{iCluster,1} = [DataForLMMFreq{iCluster,1}; CaData.(BrainStates{iState,1}).Freq(:,3)- CaData.(BrainStates{iState,1}).Freq(:,1)];
            DataForLMMAmp{iCluster,1} = [DataForLMMAmp{iCluster,1}; CaData.(BrainStates{iState,1}).Amp(:,3)- CaData.(BrainStates{iState,1}).Amp(:,1)];
            StateCond{iCluster,1} = [StateCond{iCluster,1}; repmat(iState,size(CaData.(BrainStates{iState,1}).Freq,1),1)];
            AnimalName{iCluster,1} = [AnimalName{iCluster,1}; repmat(CaData.(BrainStates{iState,1}).AnimalNameFreq,1,1)];
        end
        Cluster{iCluster,1} = [repmat(iCluster,length(AnimalName{iCluster,1}),1)];
end

%across states
T = [];
for iCluster = 1:4
    T = [T;array2table([AnimalName{iCluster,1}, DataForLMMFreq{iCluster,1}, StateCond{iCluster,1}, Cluster{iCluster,1} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Cluster = categorical(T.Cluster);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Cluster + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

T = [];
for iCluster = 1:4
    T = [T;array2table([AnimalName{iCluster,1}, DataForLMMAmp{iCluster,1}, StateCond{iCluster,1}, Cluster{iCluster,1} ],...
    'VariableNames',{'Animal','Freq', 'State', 'Cluster'})];
end
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Cluster = categorical(T.Cluster);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Cluster + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Cluster + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

hFig = figure(10023);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [100 20 500 666])
set(gcf,'color','white')
set(gcf,'name',' DiffSWSn SWSn+1','NumberTitle','off')
tiledlayout(2,2);
Labels = {'SO', 'Spindle', 'SO+Spindle','Control'};

nexttile
bar([nanmean(DiffREMTriplets.Freq{1}),nanmean(DiffREMTriplets.Freq{2}),nanmean(DiffREMTriplets.Freq{3}),nanmean(DiffREMTriplets.Freq{4})]);
[eBar] = PlotErrorbarBetween(DiffREMTriplets.Freq,PlotColor{iCluster,1});
plotSignificanceBetweenNonParametric(DiffREMTriplets.Freq)
box(gca,'off')
title({'SWS_{n}  REM  SWS_{n+1}'});
set(gca, 'XTick', 1:4,'XTickLabel', Labels);
ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca,"LineWidth",1);
set(gca, 'FontName', 'Arial')
ylim([-0.012 0])

nexttile
bar([nanmean(DiffWAKTriplets.Freq{1}),nanmean(DiffWAKTriplets.Freq{2}),nanmean(DiffWAKTriplets.Freq{3}),nanmean(DiffWAKTriplets.Freq{4})]);
[eBar] = PlotErrorbarBetween(DiffWAKTriplets.Freq,PlotColor{iCluster,1});
plotSignificanceBetweenNonParametric(DiffWAKTriplets.Freq)
box(gca,'off')
title({'SWS_{n}  Wake  SWS_{n+1}'});
set(gca, 'XTick', 1:4,'XTickLabel', Labels);
ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca,"LineWidth",1);
set(gca, 'FontName', 'Arial')
ylim([-0.012 0])

nexttile
bar([nanmean(DiffREMTriplets.Amp{1}),nanmean(DiffREMTriplets.Amp{2}),nanmean(DiffREMTriplets.Amp{3}),nanmean(DiffREMTriplets.Amp{4})]);
[eBar] = PlotErrorbarBetween(DiffREMTriplets.Amp,PlotColor{iCluster,1});
plotSignificanceBetweenNonParametric(DiffREMTriplets.Amp)
box(gca,'off')
title({'SWS_{n}  REM  SWS_{n+1}'});
set(gca, 'XTick', 1:4,'XTickLabel', Labels);
ylabel(strcat('Amplitude of Ca^{2+} transients (Hz)'))
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca,"LineWidth",1);
set(gca, 'FontName', 'Arial')
ylim([-0.06 0.06])

nexttile
bar([nanmean(DiffWAKTriplets.Amp{1}),nanmean(DiffWAKTriplets.Amp{2}),nanmean(DiffWAKTriplets.Amp{3}),nanmean(DiffWAKTriplets.Amp{4})]);
[eBar] = PlotErrorbarBetween(DiffWAKTriplets.Amp,PlotColor{iCluster,1});
plotSignificanceBetweenNonParametric(DiffWAKTriplets.Amp)
box(gca,'off')
title({'SWS_{n}  Wake  SWS_{n+1}'});
set(gca, 'XTick', 1:4,'XTickLabel', Labels);
ylabel(strcat('Amplitude of Ca^{2+} transients (Hz)'))
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca,"LineWidth",1);
set(gca, 'FontName', 'Arial')
ylim([-0.06 0.06])