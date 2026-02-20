% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
%% Triplet Analysis
close all
clear
UseIntervalBeforeEvent =1;
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

for iCluster = 1%:3
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
    CaDataState.SWS.AmpActCells = [];
    CaDataState.SWS.FreqActCells = [];
    CaDataState.REM.AmpActCells = [];
    CaDataState.REM.FreqActCells = [];
    CaDataState.Wake.AmpActCells = [];
    CaDataState.Wake.FreqActCells = [];
    CaDataState.SWS.AnimalName = [];
    CaDataState.REM.AnimalName = [];
    CaDataState.Wake.AnimalName = [];
    
    CaDataPerCellState.SWS.Amp = [];
    CaDataPerCellState.SWS.Freq = [];
    CaDataPerCellState.REM.Amp = [];
    CaDataPerCellState.REM.Freq = [];
    CaDataPerCellState.Wake.Amp = [];
    CaDataPerCellState.Wake.Freq = [];
    CaDataPerCellState.SWS.AnimalName = [];
    CaDataPerCellState.REM.AnimalName = [];
    CaDataPerCellState.Wake.AnimalName = [];
    
    for iFile = 1:length(FileName)
        load(strcat(DirData,FileName{iFile,1}),'Events','CaDataAll','SlstAll','FrameRate');
        
        events = ~isnan(CaDataAll);
        
        % Sliding window count
        winSize = 4;
        nCells = size(events,1);
        nTime = size(events,2);
        syncCounts = zeros(1, nTime-winSize+1);
        
        for t = 1:(nTime-winSize+1)
            syncCounts(t) = sum(sum(events(:,t:(t+winSize-1)),2) > 0); 
        end
        
        syncCounts(syncCounts == 1) = 0;
        syncCounts(syncCounts==0) = nan;
        fractionActive = syncCounts / nCells;
        CaDataAll = fractionActive;
        
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
            
            TMPCaDataThirdAmp = nan(3,size(TMPEpisodes,2));
            TMPCaDataThirdFre = nan(3,size(TMPEpisodes,2));
            
            for iEpisode = 1:size(TMPEpisodes,2)
                    for iThird = 1:3
                        TMPCaDataThirdAmp(iThird,iEpisode) = nanmean(nanmean(CaDataAll(:,round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate)),2));
                        TMPCaDataThirdFre(iThird,iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(iThird,iEpisode)*FrameRate)+1:round(TMPEpisodes(iThird+1,iEpisode)*FrameRate))),2)./((TMPEpisodes(2,iEpisode)-TMPEpisodes(1,iEpisode))));
                    end
            end
            
            TMPCaDataThirdAmp(isnan(TMPCaDataThirdAmp))= 0;
            TMPCaDataThirdFre(isnan(TMPCaDataThirdFre))=0;
            
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
            
            TMPCaDataAmp = nan(1,size(TMPEpisodes,2));
            TMPCaDataFre = nan(1,size(TMPEpisodes,2));
            TMPCaDataAmpActCells = nan(1,size(TMPEpisodes,2));
            TMPCaDataFreqActCells = nan(1,size(TMPEpisodes,2));
            TMPCaDataPerCellAmp = [];
            TMPCaDataPerCellFre = [];
            
            for iEpisode = 1:size(TMPEpisodes,2)
                    TMPCaDataAmp(1,iEpisode) = nanmean(nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate)),2));
                    TMPCaDataFre(1,iEpisode) = mean(sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate))),2)./(TMPEpisodes(4,iEpisode)-TMPEpisodes(1,iEpisode)));
                    TMPCaDataPerCellAmp = [TMPCaDataPerCellAmp; (nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate)),2))];
                    TMPCaDataPerCellFre = [TMPCaDataPerCellFre; (sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate))),2)./(TMPEpisodes(4,iEpisode)-TMPEpisodes(1,iEpisode)))];
                    
                    TMPActCellsFreq = (sum(~isnan(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate))),2)./(TMPEpisodes(4,iEpisode)-TMPEpisodes(1,iEpisode)));
                    TMPActCellsFreq = mean(TMPActCellsFreq(TMPActCellsFreq~=0));
                    TMPActCellsAmp = (nanmean(CaDataAll(:,round(TMPEpisodes(1,iEpisode)*FrameRate)+1:round(TMPEpisodes(4,iEpisode)*FrameRate)),2));
                    TMPActCellsAmp = mean(TMPActCellsAmp(~isnan(TMPActCellsAmp)));
                    
                    TMPCaDataAmpActCells(1,iEpisode) = TMPActCellsAmp;
                    TMPCaDataFreqActCells(1,iEpisode) = TMPActCellsFreq;
            end
            
            TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
            TMPCaDataFre(isnan(TMPCaDataFre)) = 0;
            
            switch iState
                case 1
                    CaDataState.SWS.Amp = [CaDataState.SWS.Amp; TMPCaDataAmp'];
                    CaDataState.SWS.Freq = [CaDataState.SWS.Freq; TMPCaDataFre'];
                    CaDataState.SWS.AmpActCells = [CaDataState.SWS.AmpActCells; TMPCaDataAmpActCells'];
                    CaDataState.SWS.FreqActCells = [CaDataState.SWS.FreqActCells; TMPCaDataFreqActCells'];
                    CaDataPerCellState.SWS.Amp = [CaDataPerCellState.SWS.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.SWS.Freq = [CaDataPerCellState.SWS.Freq; TMPCaDataPerCellFre];
                    CaDataState.SWS.AnimalName = [CaDataState.SWS.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.SWS.AnimalName = [CaDataPerCellState.SWS.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataPerCellFre,1),1)];
                case 2
                    CaDataState.REM.Amp = [CaDataState.REM.Amp; TMPCaDataAmp'];
                    CaDataState.REM.Freq = [CaDataState.REM.Freq; TMPCaDataFre'];
                    CaDataState.REM.AmpActCells = [CaDataState.REM.AmpActCells; TMPCaDataAmpActCells'];
                    CaDataState.REM.FreqActCells = [CaDataState.REM.FreqActCells; TMPCaDataFreqActCells'];
                    CaDataPerCellState.REM.Amp = [CaDataPerCellState.REM.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.REM.Freq = [CaDataPerCellState.REM.Freq; TMPCaDataPerCellFre];
                    CaDataState.REM.AnimalName = [CaDataState.REM.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.REM.AnimalName = [CaDataPerCellState.REM.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataPerCellFre,1),1)];
                case 3
                    CaDataState.Wake.Amp = [CaDataState.Wake.Amp; TMPCaDataAmp'];
                    CaDataState.Wake.Freq = [CaDataState.Wake.Freq; TMPCaDataFre'];
                    CaDataState.Wake.AmpActCells = [CaDataState.Wake.AmpActCells; TMPCaDataAmpActCells'];
                    CaDataState.Wake.FreqActCells = [CaDataState.Wake.FreqActCells; TMPCaDataFreqActCells'];
                    CaDataPerCellState.Wake.Amp = [CaDataPerCellState.Wake.Amp; TMPCaDataPerCellAmp];
                    CaDataPerCellState.Wake.Freq = [CaDataPerCellState.Wake.Freq; TMPCaDataPerCellFre];
                    CaDataState.Wake.AnimalName = [CaDataState.Wake.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataFre,2),1)];
                    CaDataPerCellState.Wake.AnimalName = [CaDataPerCellState.Wake.AnimalName; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataPerCellFre,1),1)];
            end
        end
    end
    
    AnimalList = [CaDataState.SWS.AnimalName; CaDataState.REM.AnimalName; CaDataState.Wake.AnimalName ];
    StateCondition = [ones(size(CaDataState.SWS.Freq));...
        ones(size(CaDataState.REM.Freq))*2;...
        ones(size(CaDataState.Wake.Freq))*3];
    CellCluster = [ones(size(CaDataState.SWS.Freq))*iCluster;...
        ones(size(CaDataState.REM.Freq))*iCluster;...
        ones(size(CaDataState.Wake.Freq))*iCluster];
    DataForAnova = [CaDataState.SWS.Amp;...
        CaDataState.REM.Amp;...
        CaDataState.Wake.Amp];
    DataForFriedman = [CaDataState.SWS.Freq, CaDataState.REM.Freq, CaDataState.Wake.Freq];
    [pValue, tbl, stats] = friedman(DataForFriedman, 1, 'off');
    
    DataForFriedman = [CaDataState.SWS.Amp, CaDataState.REM.Amp, CaDataState.Wake.Amp];
    [pValue, tbl, stats] = friedman(DataForFriedman, 1, 'off');
    
    T = [array2table([AnimalList, DataForAnova,StateCondition, CellCluster],...
        'VariableNames',{'Animal','CaDataState', 'StateCondition', 'CellClusterCondition'})];
    T.StateCondition = categorical(T.StateCondition);
    T.CellClusterCondition = categorical(T.CellClusterCondition);
    
    lm1Model = fitlme(T,'CaDataState ~ StateCondition *CellClusterCondition + (1|Animal)');
    lm2Model = fitlme(T,'CaDataState ~ StateCondition +CellClusterCondition + (1|Animal)');
    results = compare(lm2Model, lm1Model); 
    
    if results.pValue(end) > 0.05
        lm1Model = lm2Model;
    else
        PValues = results.pValue(end);
    end
    
    lm2Model = fitlme(T,'CaDataState ~ StateCondition + (1|Animal)');
    results = compare(lm2Model, lm1Model); 
    if results.pValue(end) <= 0.05
        PValuesState = results.pValue(end);
    end
    
    lm2Model = fitlme(T,'CaDataState ~ CellClusterCondition + (1|Animal)');
    results = compare(lm2Model, lm1Model); 
    if results.pValue(end) <= 0.05
        PValuesCluster = results.pValue(end);
    end
    
    %% state histograms
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
    PlotColors{1,1} = [0.8500 0.3250 0.0980];
    PlotColors{2,1} = [0 0.4470 0.7410];
    PlotColors{3,1} = [0.3 0.7 0.1];
    
    for iState = 1:3
        Data = CaDataState.(BrainStates{iState,1}).Amp*100;
        [h, edges] = histcounts(Data,[0:0.1:10]);
        h = h/sum(h)*100;
        bar(edges(2:end),h,'FaceColor',PlotColors{iState,1},'EdgeColor',[1 1 1],'FaceAlpha',0.2)
        ft = fittype( 'smoothingspline' );
        fitresult = fit(edges(2:end)', h', ft);
        hold on
        plt = plot(fitresult);
        plt.Color= PlotColors{iState,1};
        plt.LineWidth = 1.5;
    end
    
    ylim([0 20])
    xlim([0 5])
    ylabel('Probability (%)')
    xlabel('percent coactive cells')
    legend('Wake','','SWS','','REM')
    legend boxoff
    box(gca,'off')
    
    hFig = figure(33);
    if iCluster == 1
        tiledlayout(6,1)
        set(gcf,'PaperPositionMode','auto')
        set(hFig, 'Position', [0 20 500 500])
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
        Data(:,iState) = CaDataState.(BrainStates{iState,1}).Amp*100;
    end
    Data(sum(isnan(Data),2)==3,:) = [];
    Data(isnan(Data))=0;
    
    [eBar,MeanLine] = PlotErrorbar(Data,PlotColors{iCluster,1});
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca,"LineWidth",1);
    set(gca, 'FontName', 'Arial')
    ylabel(strcat('Frequency of Ca^{2+} transients (Hz)'))
    ylim([1 2])
    set(gca, 'XTick', 1:3,'XTickLabel', Labels);
    box(gca,'off')
    
    [pValue, tbl, stats] = friedman(Data, 1, 'off');
    plotSignificanceWithinNonParametric(Data)
    xlim([0.5 3.5])
    
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
        ylabel(strcat('coactive cells (%)'))
        ylim([1 2])
        xlim([0.5 3.5])
        xlabel('Third')
        plotSignificanceWithinNonParametric(Data)
    end
    box(gca,'off')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
end
 
StateCond = [];
Third = [];
DataForLMMFreq = [];
DataForLMMAmp = [];
AnimalName = [];

for iState = 1:3
    DataForLMMFreq = [DataForLMMFreq; [CaDataThirds.(BrainStates{iState,1}).Freq(:,1);...
        CaDataThirds.(BrainStates{iState,1}).Freq(:,2);...
        CaDataThirds.(BrainStates{iState,1}).Freq(:,3)]];
        
    DataForLMMAmp = [DataForLMMAmp; CaDataThirds.(BrainStates{iState,1}).Amp(:,1);...
        CaDataThirds.(BrainStates{iState,1}).Amp(:,2);...
        CaDataThirds.(BrainStates{iState,1}).Amp(:,3)];
        
    StateCond = [StateCond; repmat(iState,size(CaDataThirds.(BrainStates{iState,1}).Freq,1)*3,1)];
    
    Third = [Third; ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1);...
        ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1) *2;...
        ones(size(CaDataThirds.(BrainStates{iState,1}).Freq,1),1)*3];
        
    AnimalName = [AnimalName; repmat(CaDataThirds.(BrainStates{iState,1}).AnimalName,3,1)];
end

T = [array2table([AnimalName, DataForLMMFreq, StateCond, Third],...
    'VariableNames',{'Animal','Freq', 'State', 'Thirds'})];
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);

lm1Model = fitlme(T,strcat('Freq',' ~ State *Thirds + (1|Animal)'));
lm2Model = fitlme(T,strcat('Freq',' ~ State +Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) > 0.05
    lm1Model = lm2Model;
else
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

T = [array2table([AnimalName, DataForLMMAmp, StateCond, Third],...
    'VariableNames',{'Animal','Amp', 'State', 'Thirds'})];
T.Animal = categorical(T.Animal);
T.State = categorical(T.State);
T.Thirds = categorical(T.Thirds);

lm1Model = fitlme(T,strcat('Amp',' ~ State *Thirds + (1|Animal)'));
lm2Model = fitlme(T,strcat('Amp',' ~ State +Thirds + (1|Animal)'));
results = compare(lm2Model, lm1Model); 
if results.pValue(end) > 0.05
    lm1Model = lm2Model;
else
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