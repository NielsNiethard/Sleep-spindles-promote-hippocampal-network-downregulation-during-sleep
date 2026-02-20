% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
function [ActiveCells,InactiveCells,NumberOfEventsPerEpisode,NoEventNREMEpisodes,EpisodeDurations,EventDensity] = DefineActInactCells(CaData,FrameRate,Events,NREMEpisodes)
    
    CaDataNoEventAct = CaData;
    
    for iEvent= 1: size(Events,2)
        % set all values that are inside an event to -100 to remove later
        CaDataNoEventAct(:,round(Events(1,iEvent)*FrameRate):round(Events(2,iEvent)*FrameRate))=-100; 
    end
    
    NumberOfEventsPerEpisode = nan(1,size(NREMEpisodes,2));
    NoEventNREMEpisodes = [];
    
    for iEpisode = 1: size(NREMEpisodes,2)
        EventsPerEpisode = find(Events(1,:) >= NREMEpisodes(1,iEpisode) & Events(2,:) <= NREMEpisodes(2,iEpisode));
        CurrEventDataEvents = [];
        CurrDataOutsideEvents = CaDataNoEventAct(:,round(NREMEpisodes(1,iEpisode)*FrameRate):round(NREMEpisodes(2,iEpisode)*FrameRate));
        CurrDataOutsideEvents(:,CurrDataOutsideEvents(1,:)==-100)=[]; % remove frames during events
        
        for iEvent = 1:length(EventsPerEpisode)
            CurrEventDataEvents = [CurrEventDataEvents,...
                CaData(:,round(Events(1,EventsPerEpisode(iEvent))*FrameRate):...
                round(Events(2,EventsPerEpisode(iEvent))*FrameRate))];
        end
        
        EpisodeDurations(1,iEpisode) = NREMEpisodes(2,iEpisode) - NREMEpisodes(1,iEpisode);
        
        if ~isempty(EventsPerEpisode)
            NumberOfEventsPerEpisode(1,iEpisode) = length(EventsPerEpisode);
            EventDensity(1,iEpisode) = length(EventsPerEpisode)/EpisodeDurations(1,iEpisode);
            
            ActiveCells{1,iEpisode} = find((nansum(CurrEventDataEvents,2)./size(CurrEventDataEvents,2) -...
                nansum(CurrDataOutsideEvents,2)./size(CurrDataOutsideEvents,2)) > 0);
            
            cutoff = prctile(nansum(CurrDataOutsideEvents,2)./size(CurrDataOutsideEvents,2),20);
            InactiveCells{1,iEpisode} = find((nansum(CurrDataOutsideEvents,2)./size(CurrDataOutsideEvents,2))>0);
        else
            NoEventNREMEpisodes = [NoEventNREMEpisodes; iEpisode];
            ActiveCells{1,iEpisode} = [];
            InactiveCells{1,iEpisode} = [];
            NumberOfEventsPerEpisode(1,iEpisode) = 0;
            EventDensity(1,iEpisode) = 0;
        end
        
        InactiveCells{1,iEpisode} = setdiff(InactiveCells{1,iEpisode},ActiveCells{1,iEpisode});
    end
end