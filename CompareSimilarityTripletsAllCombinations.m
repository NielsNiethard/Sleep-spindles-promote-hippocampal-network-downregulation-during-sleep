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
        [Triplets] = ExtractSleepTripletsAllCombinations(Episodes, cfg);

        switch iCluster
            case 1
                TripletType ='WakeNREMWake';
            case 2
                TripletType ='REMWakeREM';
            case 3
                TripletType ='NREMWakeNREM';
        end
        
        Before = Triplets.(TripletType).Before;
        Middle = Triplets.(TripletType).Middle;
        After = Triplets.(TripletType).After;

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
        
        [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,Before);
        
        %%
        TMPCaDataAmp = nan(3,size(Triplets.(TripletType).Before,2));
        TMPCaDataFre = nan(3,size(Triplets.(TripletType).Before,2));
        REMEpisodeDistance = nan(1,size(Triplets.(TripletType).Before,2));
        REMEpisodeNREMNDuration = nan(1,size(Triplets.(TripletType).Before,2));
        REMEpisodeNREMN1Duration = nan(1,size(Triplets.(TripletType).Before,2));
        REMEventDensity = nan(1,size(Triplets.(TripletType).Before,2));
        
        for iTriplet = 1:3
            for iEpisode = 1:size(Triplets.(TripletType).Before,2)
                switch iTriplet
                    case 1
                        TMPEpisodes = Before;
                    case 2
                        TMPEpisodes = Middle;
                    case 3
                        TMPEpisodes = After;
                end
                
                ControlActiveCells = find(~isnan(nanmean(CaDataAll(:,round(Before(1,iEpisode)*FrameRate)+1:round(Before(2,iEpisode)*FrameRate)),2)));
                ActiveCells = find(~isnan(nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2)));
                common = intersect(ControlActiveCells, ActiveCells);
                similarity = length(common) / max(length(ControlActiveCells), length(ActiveCells)) * 100;
                
                TMPCaDataAmp(iTriplet,iEpisode) = similarity;
                TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                
                if iTriplet ==2
                    REMEpisodeDistance(1,iEpisode) = After(2,iEpisode) - Before(1,iEpisode);
                    REMEpisodeNREMNDuration(1,iEpisode) = Before(2,iEpisode) - Before(1,iEpisode);
                    REMEpisodeNREMN1Duration(1,iEpisode) = After(2,iEpisode) - After(1,iEpisode);
                    REMEventDensity(1,iEpisode) = EventDensity(1,iEpisode);
                end
            end
        end
        
        REMEventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        NumberOfEventsPerEpisode(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataFre(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
        
        EpisodeNREMNDuration.REM = [EpisodeNREMNDuration.REM , REMEpisodeNREMNDuration];
        EpisodeNREMN1Duration.REM = [EpisodeNREMN1Duration.REM, REMEpisodeNREMN1Duration];
        EpisodeDistance.REM = [EpisodeDistance.REM,REMEpisodeDistance];
        EventsPerEpisode.REM = [EventsPerEpisode.REM; NumberOfEventsPerEpisode'];
        EventsDensityPerEpisode.REM = [EventsDensityPerEpisode.REM; REMEventDensity'];
                    
        CaData.SWSREMSWS.Amp = [CaData.SWSREMSWS.Amp; TMPCaDataAmp'];
        CaData.SWSREMSWS.Freq = [CaData.SWSREMSWS.Freq; TMPCaDataFre'];
        CaData.SWSREMSWS.AnimalNameFreq = [CaData.SWSREMSWS.AnimalNameFreq; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
        CaData.SWSREMSWS.AnimalNameAmp = [CaData.SWSREMSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataAmp,2),1)];
       
        cfg = [];
        cfg.MinEpisodeDuration = MinEpisodeDuration;
        cfg.MicroarousalMerge = true;
        cfg.Mode = 'WAKE';
        
        [Triplets] = ExtractSleepTripletsAllCombinations(Episodes, cfg);

        switch iCluster
            case 1
                TripletType ='WakeREMWake';
            case 2
                TripletType ='REMNREMREM';
            case 3
                TripletType ='NREMREMNREM';
        end
        
        Before = Triplets.(TripletType).Before;
        Middle = Triplets.(TripletType).Middle;
        After = Triplets.(TripletType).After;
        
        %%
        TMPCaDataThirdAmp = nan(9,size(Triplets.(TripletType).Before,2));
        TMPCaDataThirdFre = nan(9,size(Triplets.(TripletType).Before,2));
        TMPCaDataAmp = nan(3,size(Triplets.(TripletType).Before,2));
        TMPCaDataFre = nan(3,size(Triplets.(TripletType).Before,2));
        WAKEpisodeDistance = nan(1,size(Triplets.(TripletType).Before,2));
        WAKEpisodeNREMNDuration = nan(1,size(Triplets.(TripletType).Before,2));
        WAKEpisodeNREMN1Duration = nan(1,size(Triplets.(TripletType).Before,2));
        WakeEventDensity = nan(1,size(Triplets.(TripletType).Before,2));
        
        for iTriplet = 1:3
            for iEpisode = 1:size(Triplets.(TripletType).Before,2)
                switch iTriplet
                    case 1
                        TMPEpisodes = Before;
                    case 2
                        TMPEpisodes = Middle;
                    case 3
                        TMPEpisodes = After;
                end
                
                ControlActiveCells = find(~isnan(nanmean(CaDataAll(:,round(Before(1,iEpisode)*FrameRate)+1:round(Before(2,iEpisode)*FrameRate)),2)));
                ActiveCells = find(~isnan(nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2)));
                common = intersect(ControlActiveCells, ActiveCells);
                similarity = length(common) / max(length(ControlActiveCells), length(ActiveCells)) * 100;
                
                TMPCaDataAmp(iTriplet,iEpisode) = similarity;
                TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                
                if iTriplet ==2
                    WAKEpisodeDistance(1,iEpisode) = TMPEpisodes(2,iEpisode) - TMPEpisodes(1,iEpisode);
                    WakeEventDensity(1,iEpisode) = EventDensity(1,iEpisode);
                    WAKEpisodeNREMNDuration(1,iEpisode) = Before(2,iEpisode) - Before(1,iEpisode);
                    WAKEpisodeNREMN1Duration(1,iEpisode) = After(2,iEpisode) - After(1,iEpisode);
                end
            end
        end
        
        WakeEventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        NumberOfEventsPerEpisode(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
        
        %match WAk and REM episode durations
        D=abs(WAKEpisodeDistance(:)-REMEpisodeDistance(:).');
        MatchedWakEpisodes=sortrows(matchpairs(D,max(D(:))),1);
        
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
        title(TripletType)
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
        plotSignificanceWithinNonParametric(Data)
        
        %%
        Data = CaData.SWSREMSWS.Amp;
        Data(sum((Data),2)==0,:) =[];
        Data(isnan(Data)) = 0;
        title(TripletType)
        nexttile(20*(iCluster-1)+11)
        [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
        Labels = {'SWS_{n}', 'REM', 'SWS_{n+1}'};
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 100])
        box(gca,'off')
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
        [pValue, tbl, stats] = friedman(Data, 1, 'off');
        if pValue < 0.05
            plotSignificanceWithinNonParametric(Data)
        end
        
        Data = CaData.SWSWAKSWS.Amp;
        Data(sum((Data),2)==0,:) =[];
        Data(isnan(Data)) = 0;
        title(TripletType)
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
        ylim([0 100])
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
        
        %compare SWSn+1 - SWSn betwen WAKE and REM
        title(TripletType)
        nexttile
        bar([nanmean(CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1)),nanmean(CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1))]);
        hold on
        clear Data
        Data{1,1}= CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1);
        Data{1,2}= CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1);
        [eBar] = PlotErrorbarBetween(Data,PlotColor{iCluster,1});
        Labels = {'REM', 'WAKE'};
        
        title({'SWS_{n+1}  - SWS_{n}'});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (Hz)'))
        set(gca, 'XTick', 1:2,'XTickLabel', Labels);
        ylim([-50 0.0])
        box(gca,'off')
        plotSignificanceBetweenNonParametric(Data)
        clear Data
        
        title(TripletType)
        Labels = {'REM', 'WAKE'};
        title({'SWS_{n+1}  - SWS_{n}'});
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('Amplitude of Ca^{2+} transients (Hz)'))
        set(gca, 'XTick', 1:2,'XTickLabel', Labels);
        ylim([-50 0.0])
        box(gca,'off')
        plotSignificanceBetweenNonParametric(Data)
        clear Data
        
        %correlation SWS REM SWSN+1
        CorrType = 'Spearman';
        Color = 'k';
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeDistance.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Distance SWS_{n} _{REM} SWS_{n+1}'));
        ylabel ('SWS_{n+1} - SWS_{n}')
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeDistance.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Distance SWS_{n} _{REM} SWS_{n+1}'));
        title ('Amp')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n} '));
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n} '));
        title ('Amp')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n+1}'));
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n+1}'));
        title ('Amp')
        title(TripletType)
        nexttile
        
        %correlation SWS WAK SWSN+1
        CorrType = 'Spearman';
        Color = 'k';
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeDistance.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Distance SWS_{n} _{WAK} SWS_{n+1}'));
        ylabel ('SWS_{n+1} - SWS_{n}')
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeDistance.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Distance SWS_{n} _{WAK} SWS_{n+1}'));
        title ('Amp')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n} '));
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n} '));
        title ('Amp')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n+1}'));
        title ('Freq')
        title(TripletType)
        nexttile
        [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
        xlabel(strcat('Duration SWS_{n+1}'));
        title ('Amp')
    end
    
    title(TripletType)
    nexttile
    %correlations
    CorrType = 'Spearman';
    Color = 'k';
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensityPerEpisode.REM,CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
    ylabel ('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    title(TripletType)
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensityPerEpisode.REM,CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
        ylabel ('SWS_{n+1} - SWS_{n}')
    title ('Amp')
    title(TripletType)
    nexttile
    %correlations
    CorrType = 'Spearman';
    Color = 'k';
    title(TripletType)
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensityPerEpisode.WAK,CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
        ylabel ('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    title(TripletType)
    nexttile
    [rho,pval] = PLOTCorrelations(EventsDensityPerEpisode.WAK,CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
        ylabel ('SWS_{n+1} - SWS_{n}')
    title ('Amp')
    %%
    hFig = figure(1);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
    Labels = {'1^{st}', '2^{nd}', '3^{rd}','1^{st}', '2^{nd}', '3^{rd}','1^{st}', '2^{nd}', '3^{rd}'};
    PlotColors{1,1} = [0 0.4470 0.7410];
    PlotColors{2,1} = [0.3 0.7 0.1];
    PlotColors{3,1} = [0 0.4470 0.7410];
  
    Data = CaDataThirds.SWSREMSWS.Amp;
    Data(nansum((Data),2)==0,:) = [];
    Data(isnan(Data))=0;
    title(TripletType)
    nexttile
    for iTriplet = 1:3
        for iThird = 1:3
            b = boxchart(repmat(iThird+((iTriplet-1)*3),size(Data,1),1),Data(:,iThird+((iTriplet-1)*3)),'MarkerStyle','none');
            hold all
            b.BoxWidth = 0.75;
            b.BoxFaceColor = PlotColors{iTriplet,1};
            b.BoxFaceAlpha = 0.4;
            b.LineWidth = 1.5;
            b.Notch = 'on';
        end
    end
    title({'SWS_{n}                    REM                   SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca, 'XTick', 1:9,'XTickLabel', Labels);
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    ylim([0 160])
    xlabel('Third')
    box(gca,'off')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    title(TripletType)
    
    %Mean across episodes
    hFig = figure(2);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
   
    Data = CaData.SWSREMSWS.Amp;
    Data(nansum((Data),2)==0,:) = [];
    Data(isnan(Data))=0;
    title(TripletType)
    nexttile
    for iTriplet = 1:3
        b = boxchart(repmat(iTriplet,size(Data,1),1),Data(:,iTriplet),'MarkerStyle','none');
        hold all
        b.BoxWidth = 0.75;
        b.BoxFaceColor = PlotColors{iTriplet,1};
        b.BoxFaceAlpha = 0.4;
        b.LineWidth = 1.5;
        b.Notch = 'on';
    end
    title({'SWS_{n}  REM  SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    ylim([0 100])
    xlabel('Third')
    box(gca,'off')
    title(TripletType)
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    
    %%
    %SWS WAKE SWS Triplets
    hFig = figure(3);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
    Labels = {'1^{st}', '2^{nd}', '3^{rd}','1^{st}', '2^{nd}', '3^{rd}','1^{st}', '2^{nd}', '3^{rd}'};
    PlotColors{1,1} = [0 0.4470 0.7410];
    PlotColors{2,1} = [0.8500 0.3250 0.0980];
    PlotColors{3,1} = [0 0.4470 0.7410];
    
    Data = CaDataThirds.SWSWAKSWS.Amp;
    Data(nansum((Data),2)==0,:) = [];
    Data(isnan(Data))=0;
    nexttile
    for iTriplet = 1:3
        for iThird = 1:3
            b = boxchart(repmat(iThird+((iTriplet-1)*3),size(Data,1),1),Data(:,iThird+((iTriplet-1)*3)),'MarkerStyle','none');
            hold all
            b.BoxWidth = 0.75;
            b.BoxFaceColor = PlotColors{iTriplet,1};
            b.BoxFaceAlpha = 0.4;
            b.LineWidth = 1.5;
            b.Notch = 'on';
        end
    end
    title({'SWS_{n}                    WAKE                   SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca, 'XTick', 1:9,'XTickLabel', Labels);
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    ylim([0 100])
    xlabel('Third')
    box(gca,'off')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    title(TripletType)
    
    %Mean across episodes
    hFig = figure(4);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
    Labels = {'SWS_{n}', 'WAKE', 'SWS_{n+1}'};
    PlotColors{1,1} = [0 0.4470 0.7410];
    PlotColors{2,1} = [0.8500 0.3250 0.0980];
    PlotColors{3,1} = [0 0.4470 0.7410];
    Data = CaData.SWSWAKSWS.Freq;
    Data(sum(isnan(Data),2)==3,:) = [];
    Data(isnan(Data))=0;
    nexttile
    for iTriplet = 1:3
        b = boxchart(repmat(iTriplet,size(Data,1),1),Data(:,iTriplet),'MarkerStyle','none');
        hold all
        b.BoxWidth = 0.75;
        b.BoxFaceColor = PlotColors{iTriplet,1};
        b.BoxFaceAlpha = 0.4;
        b.LineWidth = 1.5;
        b.Notch = 'on';
    end
    title({'SWS_{n}  WAKE  SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    ylim([0 0.06])
    box(gca,'off')
    title(TripletType)
    
    Data = CaData.SWSWAKSWS.Amp*100;
    Data(sum(isnan(Data),2)==3,:) = [];
    Data(isnan(Data))=0;
    nexttile
    for iTriplet = 1:3
        b = boxchart(repmat(iTriplet,size(Data,1),1),Data(:,iTriplet),'MarkerStyle','none');
        hold all
        b.BoxWidth = 0.75;
        b.BoxFaceColor = PlotColors{iTriplet,1};
        b.BoxFaceAlpha = 0.4;
        b.LineWidth = 1.5;
        b.Notch = 'on';
    end
    title({'SWS_{n}  WAKE  SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    title(TripletType)
    ylim([0 160])
    xlabel('Third')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    
    %%
    %Mean across episodes
    hFig = figure(20);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [100 20 200 1000])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
    Labels = {'SWS_{n}', 'REM', 'SWS_{n+1}'};
    PlotColors{1,1} = [0 0.4470 0.7410];
    PlotColors{2,1} = [0.3 0.7 0.1];
    PlotColors{3,1} = [0 0.4470 0.7410];
    
    Data = CaData.SWSREMSWS.Amp;
    Data(sum((Data),2)==0,:) =[];
    Data(isnan(Data)) = 0;
    nexttile
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
    title({'SWS_{n}  REM  SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    title(TripletType)
    ylim([0 100])
    box(gca,'off')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    plotSignificanceWithinNonParametric(Data)
    
    %%
    %SWS WAKE SWS Triplets
    %Mean across episodes
    hFig = figure(40);
    if iCluster ==1
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [400 20 200 1000])
        set(gcf,'color','white')
        set(gcf,'name',' EpisodeThirds','NumberTitle','off')
        tiledlayout(6,1)
    end
    Labels = {'SWS_{n}', 'WAKE', 'SWS_{n+1}'};
    PlotColors{1,1} = [0 0.4470 0.7410];
    PlotColors{2,1} = [0.8500 0.3250 0.0980];
    PlotColors{3,1} = [0 0.4470 0.7410];
    
    Data = CaData.SWSWAKSWS.Amp;
    Data(sum((Data),2)==0,:) =[];
    Data(isnan(Data)) = 0;
    nexttile
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
    title({'SWS_{n}  WAKE  SWS_{n+1}'});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    title(TripletType)
    ylim([0 100])
    xlabel('Third')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    plotSignificanceWithinNonParametric(Data)
end

BrainStates{1,1} = 'SWSREMSWS';
BrainStates{2,1} = 'SWSWAKSWS';
 
StateCond = [];
Third = [];
DataForLMMFreq = [];
DataForLMMAmp = [];
AnimalName = [];

for iState = 1:2
    DataForLMMFreq = [DataForLMMFreq; [CaData.(BrainStates{iState,1}).Freq(:,1);...
        CaData.(BrainStates{iState,1}).Freq(:,2);...
        CaData.(BrainStates{iState,1}).Freq(:,3)]];
    DataForLMMAmp = [DataForLMMAmp; CaData.(BrainStates{iState,1}).Amp(:,1);...
        CaData.(BrainStates{iState,1}).Amp(:,2);...
        CaData.(BrainStates{iState,1}).Amp(:,3)];
    StateCond = [StateCond; repmat(iState,size(CaData.(BrainStates{iState,1}).Freq,1)*3,1)];
    Third = [Third; ones(size(CaData.(BrainStates{iState,1}).Freq,1),1);...
        ones(size(CaData.(BrainStates{iState,1}).Freq,1),1) *2;...
        ones(size(CaData.(BrainStates{iState,1}).Freq,1),1)*3];
    AnimalName = [AnimalName; repmat(CaData.(BrainStates{iState,1}).AnimalNameFreq,3,1)];
end

T = [array2table([AnimalName, DataForLMMFreq, StateCond, Third],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds'})];
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Thirds + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Freq',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end

StateCond = [];
Third = [];
DataForLMMFreq = [];
DataForLMMAmp = [];
AnimalName = [];

for iState = 1:2
    DataForLMMFreq = [DataForLMMFreq; [CaData.(BrainStates{iState,1}).Amp(:,1);...
        CaData.(BrainStates{iState,1}).Amp(:,2);...
        CaData.(BrainStates{iState,1}).Amp(:,3)]];
    DataForLMMAmp = [DataForLMMAmp; CaData.(BrainStates{iState,1}).Amp(:,1);...
        CaData.(BrainStates{iState,1}).Amp(:,2);...
        CaData.(BrainStates{iState,1}).Amp(:,3)];
    StateCond = [StateCond; repmat(iState,size(CaData.(BrainStates{iState,1}).Amp,1)*3,1)];
    Third = [Third; ones(size(CaData.(BrainStates{iState,1}).Amp,1),1);...
        ones(size(CaData.(BrainStates{iState,1}).Amp,1),1) *2;...
        ones(size(CaData.(BrainStates{iState,1}).Amp,1),1)*3];
    AnimalName = [AnimalName; repmat(CaData.(BrainStates{iState,1}).AnimalNameAmp,3,1)];
end

T = [array2table([AnimalName, DataForLMMAmp, StateCond, Third],...
    'VariableNames',{'Animal','Amp', 'State', 'Thirds'})];
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);

lm1Model = fitlme(T,strcat('Amp',' ~ State *Thirds + (1|Animal)'));
lm2Model = fitlme(T,strcat('Amp',' ~ State +Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValues(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Amp',' ~ Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesDiet(1,1) = 1 - results.pValue(end);
end

lm2Model = fitlme(T,strcat('Amp',' ~ State + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) <= 0.05
    PValuesGlucose(1,1) = 1 - results.pValue(end);
end