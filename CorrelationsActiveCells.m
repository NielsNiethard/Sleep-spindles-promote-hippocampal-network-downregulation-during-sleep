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
        
        [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaDataAll,FrameRate,CurrentEvents,NREMEpisodes);
        
        %%
        TMPCaDataThirdAmp = nan(9,size(REMThirds,2));
        TMPCaDataThirdFre = nan(9,size(REMThirds,2));
        TMPCaDataAmp = nan(3,size(REMThirds,2));
        TMPCaDataFre = nan(3,size(REMThirds,2));
        REMEpisodeDistance = nan(1,size(REMThirds,2));
        REMEpisodeNREMNDuration = nan(1,size(REMThirds,2));
        REMEpisodeNREMN1Duration = nan(1,size(REMThirds,2));
        
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
                
                if iTriplet <3
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                    end
                else
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,1)*FrameRate)+1:round(TMPEpisodes(iThird+1,1)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,1)*FrameRate)+1:round(TMPEpisodes(iThird+1,1)*FrameRate))),2)./(TMPEpisodes(2,1)-TMPEpisodes(1,1)));
                    end
                end
                
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMEpisodes;
                    case 2
                        TMPEpisodes = REMEpisodes;
                    case 3
                        TMPEpisodes = NREMEpisodesAfter;
                end
                
                if iTriplet <3
                    TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2));
                    TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                else
                    TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,1)*FrameRate)+1:round(TMPEpisodes(2,1)*FrameRate)),2));
                    TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,1)*FrameRate)+1:round(TMPEpisodes(2,1)*FrameRate))),2)./(TMPEpisodes(2,1)-TMPEpisodes(1,1)));
                end
                
                if iTriplet ==2
                    REMEpisodeDistance(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    REMEpisodeNREMNDuration(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    REMEpisodeNREMN1Duration(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodesAfter(1,iEpisode);
                end
            end
        end
        
        REMEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        EventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        REMEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataFre(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        
        EpisodeNREMNDuration.REM = [EpisodeNREMNDuration.REM , REMEpisodeNREMNDuration];
        EpisodeNREMN1Duration.REM = [EpisodeNREMN1Duration.REM, REMEpisodeNREMN1Duration];
        EpisodeDistance.REM = [EpisodeDistance.REM,REMEpisodeDistance];
        EventsPerEpisode.REM = [EventsPerEpisode.REM; EventDensity'];
        
        CaDataThirds.SWSREMSWS.Amp = [CaDataThirds.SWSREMSWS.Amp; TMPCaDataThirdAmp'];
        CaDataThirds.SWSREMSWS.Freq = [CaDataThirds.SWSREMSWS.Freq; TMPCaDataThirdFre'];
        CaData.SWSREMSWS.Amp = [CaData.SWSREMSWS.Amp; TMPCaDataAmp'];
        CaData.SWSREMSWS.Freq = [CaData.SWSREMSWS.Freq; TMPCaDataFre'];
        
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
                
                if iTriplet <3
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                    end
                else
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird+((iTriplet-1)*3),iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,1)*FrameRate)+1:round(TMPEpisodes(iThird+1,1)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird+((iTriplet-1)*3),iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(iThird,1)*FrameRate)+1:round(TMPEpisodes(iThird+1,1)*FrameRate))),2)./(TMPEpisodes(2,1)-TMPEpisodes(1,1)));
                    end
                end
                
                switch iTriplet
                    case 1
                        TMPEpisodes = NREMEpisodes;
                    case 2
                        TMPEpisodes = WAKEpisodes;
                    case 3
                        TMPEpisodes = NREMEpisodesAfter;
                end
                
                if iTriplet <3
                    TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate)),2));
                    TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(2,iEpisode)*FrameRate))),2)./(TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode)));
                else
                    TMPCaDataAmp(iTriplet,iEpisode) = nanmean(nanmean(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,1)*FrameRate)+1:round(TMPEpisodes(2,1)*FrameRate)),2));
                    TMPCaDataFre(iTriplet,iEpisode) = mean(sum(~isnan(CaDataAll(ActiveCells{1,iEpisode},round(TMPEpisodes(1,1)*FrameRate)+1:round(TMPEpisodes(2,1)*FrameRate))),2)./(TMPEpisodes(2,1)-TMPEpisodes(1,1)));
                end
                
                if iTriplet ==2
                    WAKEpisodeDistance(1,iEpisode) = TMPEpisodes(2,iEpisode) - TMPEpisodes(1,iEpisode);
                    WAKEpisodeNREMNDuration(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
                    WAKEpisodeNREMN1Duration(1,iEpisode) = NREMEpisodesAfter(2,iEpisode) - NREMEpisodesAfter(1,iEpisode);
                end
            end
        end
        
        WAKEpisodeNREMNDuration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeNREMN1Duration(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        EventDensity(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        WAKEpisodeDistance(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataThirdFre(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataThirdAmp(:,sum(isnan(TMPCaDataThirdAmp),1)==size(TMPCaDataThirdAmp,1)) = [];
        TMPCaDataFre(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        
        TMPCaDataThirdAmp(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataThirdFre(isnan(TMPCaDataThirdAmp)) = 0;
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
        
        CaDataThirds.SWSWAKSWS.Amp = [CaDataThirds.SWSWAKSWS.Amp; TMPCaDataThirdAmp'];
        CaDataThirds.SWSWAKSWS.Freq = [CaDataThirds.SWSWAKSWS.Freq; TMPCaDataThirdFre'];
        
        %match WAk and REM episode durations
        D=abs(WAKEpisodeDistance(:)-REMEpisodeDistance(:).');
        MatchedWakEpisodes=sortrows(matchpairs(D,max(D(:))),1);
        
        EpisodeNREMNDuration.WAK = [EpisodeNREMNDuration.WAK , WAKEpisodeNREMNDuration(:,MatchedWakEpisodes(:,1))];
        EpisodeNREMN1Duration.WAK = [EpisodeNREMN1Duration.WAK, WAKEpisodeNREMN1Duration(:,MatchedWakEpisodes(:,1))];
        CaData.SWSWAKSWS.Amp = [CaData.SWSWAKSWS.Amp; TMPCaDataAmp(:,MatchedWakEpisodes(:,1))'];
        CaData.SWSWAKSWS.Freq = [CaData.SWSWAKSWS.Freq; TMPCaDataFre(:,MatchedWakEpisodes(:,1))'];
        EpisodeDistance.WAK = [EpisodeDistance.WAK, WAKEpisodeDistance(1,MatchedWakEpisodes(:,1))];
        EventsPerEpisode.WAK = [EventsPerEpisode.WAK; EventDensity(1,MatchedWakEpisodes(:,1))'];
    end
    
    %%
    %Mean across episodes
    hFig = figure(iCluster);
    set(gcf,'PaperPositionMode','auto')
    set(hFig, 'Position', [100 20 1600 1000])
    set(gcf,'color','white')
    set(gcf,'name',strcat(' Summary',CellClusterName{iCluster,1},' Cells'),'NumberTitle','off')
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
    ylim([0 0.06])
    box(gca,'off')
    plotSignificanceWithinTTest(Data)
    
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
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    ylim([0 0.06])
    box(gca,'off')
    plotSignificanceWithinTTest(Data)
    
    %stats
    % Clean up Data1
    Data1 = CaData.SWSREMSWS.Freq;
    Data1(isnan(Data1)) = 0;  % Replace NaNs with 0
    Data1 = Data1(~any(isnan(Data1), 2), :);  % Remove rows that have NaNs
    % Clean up Data2
    Data2 = CaData.SWSWAKSWS.Freq;
    Data2(isnan(Data2)) = 0;  % Replace NaNs with 0
    Data2 = Data2(~any(isnan(Data2), 2), :);  % Remove rows that have NaNs
    
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
    plotSignificanceWithinTTest(Data)
    
    Data = CaData.SWSWAKSWS.Amp*100;
    Data(sum(isnan(Data),2)==9,:) =[];
    Data(isnan(Data)) = 0;
    
    nexttile(20*(iCluster-1)+12)
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColor{iCluster,1});
    Labels = {'SWS_{n}', 'WAKE', 'SWS_{n+1}'};
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Amplitude of Ca^{2+} transients (%','\DeltaF','/F)'))
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    plotSignificanceWithinTTest(Data)
    ylim([0 160])
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    
    %stats
    % Clean up Data1
    Data1 = CaData.SWSREMSWS.Amp;
    Data1(isnan(Data1)) = 0;  % Replace NaNs with 0
    Data1 = Data1(~any(isnan(Data1), 2), :);  % Remove rows that have NaNs
    % Clean up Data2
    Data2 = CaData.SWSWAKSWS.Amp;
    Data2(isnan(Data2)) = 0;  % Replace NaNs with 0
    Data2 = Data2(~any(isnan(Data2), 2), :);  % Remove rows that have NaNs
    
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
    
    %correlation SWS REM SWSN+1
    CorrType = 'Spearman';
    Color = 'k';
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeDistance.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Distance SWS_{n} _{REM} SWS_{n+1}'));
    ylabel ('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeDistance.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Distance SWS_{n} _{REM} SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n} '));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n} '));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.REM',CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.REM',CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    %correlation SWS WAK SWSN+1
    CorrType = 'Spearman';
    Color = 'k';
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeDistance.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Distance SWS_{n} _{WAK} SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeDistance.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Distance SWS_{n} _{WAK} SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n} '));
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMNDuration.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n} '));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.WAK',CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EpisodeNREMN1Duration.WAK',CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('Duration SWS_{n+1}'));
    ylabel('SWS_{n+1} - SWS_{n}');
    title ('Amp')
    
    nexttile
    %correlations
    CorrType = 'Spearman';
    Color = 'k';
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsPerEpisode.REM,CaData.SWSREMSWS.Freq(:,3)-CaData.SWSREMSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
    ylabel('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsPerEpisode.REM,CaData.SWSREMSWS.Amp(:,3)-CaData.SWSREMSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
    ylabel('SWS_{n+1} - SWS_{n}')
    title ('Amp')
    
    nexttile
    %correlations
    CorrType = 'Spearman';
    Color = 'k';
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsPerEpisode.WAK,CaData.SWSWAKSWS.Freq(:,3)-CaData.SWSWAKSWS.Freq(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
    ylabel('SWS_{n+1} - SWS_{n}')
    title ('Freq')
    
    nexttile
    [rho,pval] = PLOTCorrelations(EventsPerEpisode.WAK,CaData.SWSWAKSWS.Amp(:,3)-CaData.SWSWAKSWS.Amp(:,1),CorrType,Color);
    xlabel(strcat('# ',CellClusterName{iCluster,1}));
    ylabel('SWS_{n+1} - SWS_{n}')
    title ('Amp')
end