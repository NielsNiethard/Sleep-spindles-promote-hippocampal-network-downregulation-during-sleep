% This code was made by Niels Niethard niels.niethard@uni-tuebingen.de
function [Doublets] = ExtractSleepDoublets(Episodes, cfg)
% ExtractSleepDoublets
% Extracts same-state doublets:
%   NREM–NREM
%   REM–REM
%   WAKE–WAKE
%
% Each doublet is stored in Doublets.<State> with fields:
%   .First   - first episode of the pair
%   .Second  - second episode of the pair
%
% Example:
%   Doublets.NREM.First   -> NREM_n
%   Doublets.NREM.Second  -> NREM_{n+1}

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
% STEP 2: Build doublets for each state
% -------------------------------------------------------------------------
Doublets.NREM = build_doublets(NREMEpisodes);
Doublets.REM  = build_doublets(REMEpisodes);
Doublets.WAKE = build_doublets(WAKEpisodes);

end

% -------------------------------------------------------------------------
%% Helper Function: Build same-state doublets
% -------------------------------------------------------------------------
function D = build_doublets(StateEpisodes)
% Build doublets of type: State_n – State_{n+1}
D = struct();
State = StateEpisodes;
FirstIdx  = 1:(size(State,2)-1);
SecondIdx = 2:size(State,2);

% Ensure the second episode starts after the first one ends
valid = State(1,SecondIdx) > State(2,FirstIdx);
FirstIdx  = FirstIdx(valid);
SecondIdx = SecondIdx(valid);

D.First  = State(:, FirstIdx);
D.Second = State(:, SecondIdx);
end