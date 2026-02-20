% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
function [Triplets] = ExtractSleepTriplets(Episodes, cfg)
    % Extract WAKE–NREM–REM–NREM triplets with optional microarousal merging
    
    % Unpack config
    MinEpisodeDuration = cfg.MinEpisodeDuration;
    MicroarousalMerge = cfg.MicroarousalMerge;
    
    if isfield(cfg,'Mode')
        Mode = upper(cfg.Mode); % 'REM' or 'WAKE'
    else
        Mode ='REM';
    end
    
    MaxMicroarousalGap = 20; % seconds
    
    % Helper: Filter episodes by minimum duration
    filter_by_duration = @(ep) (ep(2,:) - ep(1,:) >= MinEpisodeDuration);
    
    % Helper: Recompute thirds
    compute_thirds = @(episodes) arrayfun(@(i) ...
        get_thirds(episodes(1,i), episodes(2,i)), 1:size(episodes,2), 'UniformOutput', false);
        
    %% STEP 1: Handle optional microarousal merging in NREM
    TMPEpisodes = Episodes.NREMEpisodes;
    if MicroarousalMerge
        EpisodesNew = TMPEpisodes;
        Interruption = diff(EpisodesNew(1,2:end) - EpisodesNew(2,1:end-1) - 1);
        MicroIdx = find(Interruption <= MaxMicroarousalGap);
        for i = 1:length(MicroIdx)
            EpisodesNew(2,MicroIdx(i)) = EpisodesNew(2,MicroIdx(i)+1);
            EpisodesNew(1,MicroIdx(i)+1) = nan;
        end
        EpisodesNew(:,isnan(EpisodesNew(1,:))) = [];
        TMPEpisodes = EpisodesNew;
    end
    
    % Recalculate thirds after potential merging
    ThirdsCell = compute_thirds(TMPEpisodes);
    TMPEpisodesThirds = cell2mat(ThirdsCell);
    
    %% STEP 2: Filter episodes by duration
    validNREM = filter_by_duration(TMPEpisodes);
    NREMEpisodes = TMPEpisodes(:,validNREM);
    NREMThirds = TMPEpisodesThirds(:,validNREM);
    
    validWAKE = filter_by_duration(Episodes.WAKEpisodes);
    WAKEpisodes = Episodes.WAKEpisodes(:,validWAKE);
    WAKThirds = Episodes.WAKEpisodesThirds(:,validWAKE);
    
    REMEpisodes = Episodes.REMEpisodes;
    REMThirds = Episodes.REMEpisodesThirds;
    
    %% STEP 3: Match each REM to preceding and following NREM
    %% Mode switch: REM-anchored or WAKE-anchored
    switch Mode
        case 'REM'
            BeforeNREM = nan(1, size(REMEpisodes,2));
            AfterNREM = nan(1, size(REMEpisodes,2));
            for i = 1:size(REMEpisodes,2)
                REMStart = REMEpisodes(1,i);
                pre = find(NREMEpisodes(2,:) < REMStart, 1, 'last');
                post = find(NREMEpisodes(1,:) > REMEpisodes(2,i), 1, 'first');
                if ~isempty(pre) && ~isempty(post)
                    BeforeNREM(i) = pre;
                    AfterNREM(i) = post;
                end
            end
            
            validREM = ~isnan(BeforeNREM);
            REMEpisodes = REMEpisodes(:,validREM);
            REMThirds = REMThirds(:,validREM);
            BeforeNREM = BeforeNREM(validREM);
            AfterNREM = AfterNREM(validREM);
            
            %% STEP 4: Find WAKE episode before the first NREM in each triplet
            WAKEbeforeNREM = nan(1, length(BeforeNREM));
            for i = 1:length(BeforeNREM)
                NREMstart = NREMEpisodes(1, BeforeNREM(i));
                idx = find(WAKEpisodes(2,:) < NREMstart, 1, 'last');
                if ~isempty(idx)
                    WAKEbeforeNREM(i) = idx;
                end
            end
            
            % Final filtering
            validWAKE = ~isnan(WAKEbeforeNREM);
            WAKEbeforeNREM = WAKEbeforeNREM(validWAKE);
            BeforeNREM = BeforeNREM(validWAKE);
            AfterNREM = AfterNREM(validWAKE);
            REMEpisodes = REMEpisodes(:,validWAKE);
            REMThirds = REMThirds(:,validWAKE);
            
            %% STEP 5: Assemble triplets
            Triplets.WAKEBef = WAKEpisodes(:, WAKEbeforeNREM);
            Triplets.WAKEBefThirds = WAKThirds(:, WAKEbeforeNREM);
            
        case 'WAKE'
            BeforeNREM = nan(1, size(WAKEpisodes,2));
            AfterNREM = nan(1, size(WAKEpisodes,2));
            for i = 1:size(WAKEpisodes,2)
                REMStart = WAKEpisodes(1,i);
                pre = find(NREMEpisodes(2,:) < REMStart, 1, 'last');
                post = find(NREMEpisodes(1,:) > WAKEpisodes(2,i), 1, 'first');
                if ~isempty(pre) && ~isempty(post)
                    BeforeNREM(i) = pre;
                    AfterNREM(i) = post;
                end
            end
            
            validWAK = ~isnan(BeforeNREM);
            WAKEpisodesInterl = WAKEpisodes(:,validWAK);
            WakThirdsInterl = WAKThirds(:,validWAK);
            BeforeNREM = BeforeNREM(validWAK);
            AfterNREM = AfterNREM(validWAK);
            
            Triplets.WAKE = WAKEpisodesInterl;
            Triplets.WAKEThirds = WakThirdsInterl;
            
        otherwise
            error('cfg.Mode must be ''REM'' or ''WAKE''');
    end
    
    Triplets.NREM_before = NREMEpisodes(:, BeforeNREM);
    Triplets.NREM_beforeThirds = NREMThirds(:, BeforeNREM);
    Triplets.REM = REMEpisodes;
    Triplets.REMThirds = REMThirds;
    Triplets.NREM_after = NREMEpisodes(:, AfterNREM);
    Triplets.NREM_afterThirds = NREMThirds(:, AfterNREM);
end

%% Helper Function: Split episode into thirds
function thirds = get_thirds(startT, endT)
    len = endT - startT + 1;
    t1 = floor(len / 3);
    t2 = floor(len / 3);
    t3 = len - t1 - t2;
    thirds = [
        startT;
        startT + t1 - 1;
        startT + t1 + t2 - 1;
        endT
    ];
end