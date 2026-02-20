function [Episodes] = SleepEpisodeDetector(cfg)

%detect brain state episodes
% NREM
NREMBegEpisode = strfind(any(cfg.scoring(:,1)==cfg.code_NREM,2)',[0 1]); % where does scoring flip to NREM
NREMEndEpisode = strfind(any(cfg.scoring(:,1)==cfg.code_NREM,2)',[1 0]); % where does scoring flip from NREM to something else
NREMBegEpisode = NREMBegEpisode+1; % because it always finds the epoch before
if any(cfg.scoring(1,1)==cfg.code_NREM,2)
    NREMBegEpisode = [1 NREMBegEpisode];
end
if any(cfg.scoring(end,1)==cfg.code_NREM,2)
    NREMEndEpisode = [NREMEndEpisode length(cfg.scoring)];
end
NREMEpisodes = [(NREMBegEpisode-1)*cfg.scoring_epoch_length+1; NREMEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec

% REM
if ~isempty(cfg.code_REM)
    REMBegEpisode = strfind(any(cfg.scoring(:,1)==cfg.code_REM,2)',[0 1]);
    REMEndEpisode = strfind(any(cfg.scoring(:,1)==cfg.code_REM,2)',[1 0]);
    REMBegEpisode = REMBegEpisode+1;
    if any(cfg.scoring(1,1)==cfg.code_REM,2)
        REMBegEpisode = [1 REMBegEpisode];
    end
    if any(cfg.scoring(end,1)==cfg.code_REM,2)
        REMEndEpisode = [REMEndEpisode length(cfg.scoring)];
    end
    REMEpisodes = [(REMBegEpisode-1)*cfg.scoring_epoch_length+1; REMEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec
else
    REMEpisodes = [];
end

% Wake
WAKBegEpisode = strfind((cfg.scoring(:,1)==cfg.code_WAKE)',[0 1]);
WAKEndEpisode = strfind((cfg.scoring(:,1)==cfg.code_WAKE)',[1 0]);
WAKBegEpisode = WAKBegEpisode+1;
if cfg.scoring(1,1) == cfg.code_WAKE
    WAKBegEpisode = [1 WAKBegEpisode];
end
if cfg.scoring(end,1) == cfg.code_WAKE
    WAKEndEpisode = [WAKEndEpisode length(cfg.scoring)];
end
WAKEpisodes = [(WAKBegEpisode-1)*cfg.scoring_epoch_length+1; WAKEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec


Episodes.NREMEpisodes = NREMEpisodes;
Episodes.REMEpisodes = REMEpisodes;
Episodes.WAKEpisodes = WAKEpisodes;


for iState = 1:3
    switch iState
        case 1
            TMPEpisodes = Episodes.NREMEpisodes;
        case 2
            TMPEpisodes = Episodes.REMEpisodes;
        case 3
            TMPEpisodes = Episodes.WAKEpisodes;
    end
    EpisodeDurations = TMPEpisodes(2,:) -TMPEpisodes (1,:)+1;
    EpisodeThirds = nan(4,size(EpisodeDurations,2));
    for iEpisode = 1:size(TMPEpisodes,2)
        episodeStart = TMPEpisodes(1,iEpisode);
        episodeEnd   = TMPEpisodes(2,iEpisode);
        episodeLen   = episodeEnd - episodeStart + 1;

        % Compute thirds with better handling of rounding
        third1Len = floor(episodeLen / 3);
        third2Len = floor(episodeLen / 3);
        third3Len = episodeLen - third1Len - third2Len;

        % Thirds endpoints
        t1_end = episodeStart + third1Len - 1;
        t2_end = t1_end + third2Len;
        t3_end = episodeEnd;

        % Save
        EpisodeThirds(:,iEpisode) = [episodeStart; t1_end; t2_end; t3_end];
    end
    switch iState
        case 1
            Episodes.NREMEpisodesThirds = EpisodeThirds;
        case 2
            Episodes.REMEpisodesThirds = EpisodeThirds;
        case 3
            Episodes.WAKEpisodesThirds = EpisodeThirds;
    end
end

end