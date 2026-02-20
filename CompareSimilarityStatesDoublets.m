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
        [Doublets] = ExtractSleepDoublets(Episodes, cfg);
     
        %%
        switch iCluster
            case 1
                TMPEpisodes = Doublets.WAKE;
            case 2
                TMPEpisodes = Doublets.NREM;
            case 3
                TMPEpisodes = Doublets.REM;
        end
        
        TMPCaDataAmp = nan(1,size(TMPEpisodes.First,2));
        
        for iTriplet = 1
            for iEpisode = 1:size(TMPEpisodes.First,2)
                ControlActiveCells = find(~isnan(nanmean(CaDataAll(:,round(TMPEpisodes.First(1,iEpisode)*FrameRate)+1:round(TMPEpisodes.First(2,iEpisode)*FrameRate)),2)));
                ActiveCells = find(~isnan(nanmean(CaDataAll(:,round(TMPEpisodes.Second(1,iEpisode)*FrameRate)+1:round(TMPEpisodes.Second(2,iEpisode)*FrameRate)),2)));
                common = intersect(ControlActiveCells, ActiveCells);
                similarity = length(common) / (length(ControlActiveCells)) * 100;
                TMPCaDataAmp(1,iEpisode) = similarity;
            end
        end
        
        TMPCaDataAmp(:,sum(isnan(TMPCaDataAmp),1)==size(TMPCaDataAmp,1)) = [];
        TMPCaDataAmp(isnan(TMPCaDataAmp)) = 0;
        
        CaData.SWSREMSWS.Amp = [CaData.SWSREMSWS.Amp; TMPCaDataAmp'];
        CaData.SWSREMSWS.AnimalNameAmp = [CaData.SWSREMSWS.AnimalNameAmp; repmat(str2double(FileName{iFile}(2:3)),size(TMPCaDataAmp,2),1)];
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
        Labels = {'Wake', 'SWS', 'REM'};
        PlotColors{1,1} = [0 0.4470 0.7410];
        PlotColors{2,1} = [0.3 0.7 0.1];
        PlotColors{3,1} = [0 0.4470 0.7410];
    end
    
    Data{1,iCluster} = CaData.SWSREMSWS.Amp;
    
    nexttile(20*(iCluster-1)+1)
    if iCluster ==3
        for i = 1:3
            bar (i,mean(Data{1,i}));
            hold all
        end
        hold all
        [eBar] = PlotErrorbarBetween(Data);
        set(gca,'FontSize',8,'TickLength',[0.025 0.025])
        set(gca,'TickDir','out');
        set(gca,"LineWidth",1);
        set(gca, 'FontName', 'Arial')
        ylabel(strcat('overlap active cells (%)'))
        set(gca, 'XTick', 1:3,'XTickLabel', Labels);
        ylim([0 90])
        box(gca,'off')
        
        group = [repmat({'Wake'}, 1, numel(Data{1})), ...
            repmat({'SWS'}, 1, numel(Data{2})), ...
            repmat({'REM'}, 1, numel(Data{3}))];
        data = [Data{1,1}',Data{1,2}',Data{1,3}'];
        [pValue, tbl, stats] = kruskalwallis(data, group);
        
        if pValue < 0.05
            plotSignificanceBetweenNonParametric(Data)
        end
        title ('Epoch_n - Epoch_n_+_1')
    end
end

BrainStates{1,1} = 'SWSREMSWS';
BrainStates{2,1} = 'SWSWAKSWS';
 
StateCond = [];
Third = [];
DataForLMMFreq = [];
DataForLMMAmp = [];
AnimalName = [];

for iState = 1
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