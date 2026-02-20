% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
function [Triplets] = ExtractSleepTripletsAllCombinations(Episodes, cfg)
    % ExtractSleepTriplets
    % Extracts cross-state triplets such as:
    %   WAKE–NREM–WAKE
    %   WAKE–REM–WAKE
    %   REM–NREM–REM
    %   REM–WAKE–REM
    %   NREM–WAKE–NREM
    %   NREM–REM–NREM
    %
    % Each triplet is stored in Triplets.<Pattern> with fields:
    %   .Before   - first state episode
    %   .Middle   - intervening episode
    %   .After    - second state episode
    
    % -------------------------------------------------------------------------
    % CONFIG
    % -------------------------------------------------------------------------
    MinEpisodeDuration = cfg.MinEpisodeDuration;
    filter_by_duration = @(ep) (ep(2,:) - ep(1,:) >= MinEpisodeDuration);
    
    % -------------------------------------------------------------------------
    % STEP 1: Filter episodes by duration
    % -------------------------------------------------------------------------
    validNREM = filter_by_duration(Episodes.NREMEpisodes);
    NREMEpisodes = Episodes.NREMEpisodes(:, validNREM);
    
    validWAKE = filter_by_duration(Episodes.WAKEpisodes);
    WAKEpisodes = Episodes.WAKEpisodes(:, validWAKE);
    
    validREM = filter_by_duration(Episodes.REMEpisodes);
    REMEpisodes = Episodes.REMEpisodes(:, validREM);
    
    % -------------------------------------------------------------------------
    % STEP 2: Build all cross-state triplets
    % -------------------------------------------------------------------------
    Triplets.WakeNREMWake = build_cross_triplets(WAKEpisodes, NREMEpisodes);
    Triplets.WakeREMWake  = build_cross_triplets(WAKEpisodes, REMEpisodes);
    Triplets.REMNREMREM   = build_cross_triplets(REMEpisodes, NREMEpisodes);
    Triplets.REMWakeREM   = build_cross_triplets(REMEpisodes, WAKEpisodes);
    Triplets.NREMWakeNREM = build_cross_triplets(NREMEpisodes, WAKEpisodes);
    Triplets.NREMREMNREM  = build_cross_triplets(NREMEpisodes, REMEpisodes);
end

% -------------------------------------------------------------------------
%% Helper Function: Build cross-state triplets
% -------------------------------------------------------------------------
function CrossTriplets = build_cross_triplets(StateEpisodes, MiddleEpisodes)
    % Build triplets of type: State–Middle–State (e.g., WAKE–NREM–WAKE)
    CrossTriplets = struct();
    State1 = StateEpisodes;
    State2 = MiddleEpisodes;
    
    BeforeIdx = nan(1, size(State1,2));
    MiddleIdx = nan(1, size(State1,2));
    AfterIdx  = nan(1, size(State1,2));
    
    for i = 1:size(State1,2)-1
        thisEnd = State1(2,i);
        nextStart = State1(1,i+1);
        
        % Find if a Middle episode occurs between
        mid = find(State2(1,:) > thisEnd & State2(2,:) < nextStart, 1, 'first');
        if ~isempty(mid)
            BeforeIdx(i) = i;
            MiddleIdx(i) = mid;
            AfterIdx(i)  = i + 1;
        end
    end
    
    valid = ~isnan(MiddleIdx);
    BeforeIdx = BeforeIdx(valid);
    MiddleIdx = MiddleIdx(valid);
    AfterIdx  = AfterIdx(valid);
    
    CrossTriplets.Before = State1(:, BeforeIdx);
    CrossTriplets.Middle = State2(:, MiddleIdx);
    CrossTriplets.After  = State1(:, AfterIdx);
end