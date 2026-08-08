function [response, diagnostics] = measurePeriPeakResponse(pupil, pupilTime, ...
    peakTimes, cfg)
%MEASUREPERIPEAKRESPONSE Discovery-identical peri-peak pupil measurement.
%
% This is a port of the discovery cohort's measurement, not a reimplementation
% of its description. Where the prose in docs/replication_plan_ebrains.md and
% the code in PupilHG.m disagree, the code wins, because the point of the
% exercise is that a difference in result must be attributable to the data
% rather than to the measurement. Two such disagreements are worth naming:
%
%   * The plan says the pupil is "baselined -20 to -10 s". The discovery code
%     does no such baseline. It linearly detrends each trial across the whole
%     +/-20 s window (PupilHG.m:646), averages, then subtracts the median of
%     the average (PupilHG.m:478). That is what happens here.
%   * The plan says the surrogate band comes from "100 random peak-time
%     shifts". The discovery code draws peak times uniformly at random from the
%     valid pupil range (PupilHG.m:630), which is not a shift of the observed
%     times. That is what happens here.
%
% The permutation structure is likewise copied: for each of NUMPERM
% permutations, up to 1000 observed peak times are resampled with replacement,
% an equal number of uniformly random times is drawn, trials whose standard
% deviation is an outlier are discarded, each trial is linearly detrended, and
% the surviving trials are averaged. The band is 1.96 times the standard
% deviation *across permutations*, and the surrogate band is summarised by its
% RMS over the whole window (PupilHG.m:506).

arguments
    pupil (:,1) double
    pupilTime (:,1) double
    peakTimes (:,1) double
    cfg (1,1) struct
end

response = struct('RespSig', 0, 'RespAreaNet', 0, 'RespAreaAbs', 0, ...
    'NPeaks', numel(peakTimes), 'NPeaksUsed', 0, 'MeanTrialCount', 0, ...
    'WindowMissingFraction', NaN, 'EventOverlapFraction', NaN);
diagnostics = struct('t', [], 'pd', [], 'pdse', [], 'pdrand', [], ...
    'pdserand', [], 'surrogateLevel', NaN);

fs = cfg.replication.pupilFs;
lagSeconds = (cfg.replication.timeRangeMs(1):1000 / fs:...
    cfg.replication.timeRangeMs(2))' / 1000;
nLag = numel(lagSeconds);
if nLag < 3 || isempty(peakTimes)
    return
end

% Uniform grid assumption: the staged pupil trace is a contiguous 150 Hz
% series, so a time can be converted to an index arithmetically.
t0 = pupilTime(1);
toIndex = @(t) round((t - t0) * fs) + 1;

% Valid range excludes the interpolated edges, matching obj.ValidPDRng.
finite = isfinite(pupil);
firstValid = find(finite, 1, 'first');
lastValid = find(finite, 1, 'last');
if isempty(firstValid)
    return
end

lagOffsets = round(lagSeconds * fs);
peakIndex = toIndex(peakTimes);
keep = (peakIndex + lagOffsets(1) >= firstValid) & ...
    (peakIndex + lagOffsets(end) <= lastValid);
peakIndex = peakIndex(keep);
if isempty(peakIndex)
    return
end
response.NPeaksUsed = numel(peakIndex);

% Reject windows that are mostly missing pupil, prespecified at 40%.
windowRows = peakIndex' + lagOffsets;
missingByWindow = mean(~finite(windowRows), 1)';
usable = missingByWindow <= cfg.replication.maxWindowMissingFraction;
response.WindowMissingFraction = mean(missingByWindow);
peakIndex = peakIndex(usable);
if numel(peakIndex) < cfg.replication.minPeaksPerContact
    response.NPeaksUsed = numel(peakIndex);
    return
end
response.NPeaksUsed = numel(peakIndex);

nPerm = cfg.replication.numPermutations;
maxTrials = cfg.replication.maxTrialsPerPermutation;
observedMean = nan(nLag, nPerm);
surrogateMean = nan(nLag, nPerm);
trialCount = nan(1, nPerm);

lowestStart = firstValid - lagOffsets(1);
highestStart = lastValid - lagOffsets(end);

filled = pupil;
filled(~finite) = NaN;

for k = 1:nPerm
    nDraw = min(numel(peakIndex), maxTrials);
    drawn = peakIndex(randi(numel(peakIndex), [nDraw, 1]));
    randomStart = randi([lowestStart, highestStart], [nDraw, 1]);

    observedMean(:, k) = localPermutationMean(filled, drawn, lagOffsets);
    surrogateMean(:, k) = localPermutationMean(filled, randomStart, lagOffsets);
    trialCount(k) = nDraw;
end

pd = mean(observedMean, 2, 'omitnan');
pd = pd - median(pd, 'omitnan');
pdse = std(observedMean, 0, 2, 'omitnan') * 1.96;
pdrand = mean(surrogateMean, 2, 'omitnan');
pdrand = pdrand - median(pdrand, 'omitnan');
pdserand = std(surrogateMean, 0, 2, 'omitnan') * 1.96;

if all(isnan(pd)) || all(isnan(pdse))
    return
end

% Feature extraction is restricted to +/-5 s, but the surrogate level is the
% RMS over the entire window, exactly as in the discovery code.
inWindow = lagSeconds > -cfg.replication.responseWindowSeconds & ...
    lagSeconds < cfg.replication.responseWindowSeconds;
surrogateLevel = rms(pdserand(isfinite(pdserand)));

upper = pd(inWindow) + pdse(inWindow);
lower = pd(inWindow) - pdse(inWindow);

[posFraction, posArea] = localLongestRun(lower > surrogateLevel, ...
    lower, surrogateLevel, +1);
[negFraction, negArea] = localLongestRun(upper < -surrogateLevel, ...
    upper, surrogateLevel, -1);

response.RespAreaNet = sum([posArea; negArea]);
response.RespAreaAbs = sum(abs([posArea; negArea]));
response.RespSig = max(posFraction, negFraction);
response.MeanTrialCount = mean(trialCount);

diagnostics.t = lagSeconds;
diagnostics.pd = pd;
diagnostics.pdse = pdse;
diagnostics.pdrand = pdrand;
diagnostics.pdserand = pdserand;
diagnostics.surrogateLevel = surrogateLevel;
end

% -------------------------------------------------------------------------
function averaged = localPermutationMean(pupil, starts, lagOffsets)
%LOCALPERMUTATIONMEAN Trial matrix, outlier rejection, detrend, average.

rows = starts' + lagOffsets;
trials = pupil(rows);

% Discard trials whose variability is an outlier, then remove a linear trend
% from each. detrend cannot take NaN, so gaps are filled per trial first.
spread = std(trials, 0, 1, 'omitnan');
trials(:, isoutlier(spread)) = [];
if isempty(trials)
    averaged = nan(numel(lagOffsets), 1);
    return
end

% fillmissing and detrend both operate column-wise on a matrix, so the whole
% trial set is conditioned in two calls rather than in a loop over up to a
% thousand columns. Trials that are entirely missing are restored to NaN
% afterwards, since fillmissing leaves them NaN but detrend would not.
missing = ~isfinite(trials);
if any(missing, 'all')
    entirelyMissing = all(missing, 1);
    trials = fillmissing(trials, 'linear', 'EndValues', 'nearest');
    trials = detrend(trials);
    trials(:, entirelyMissing) = NaN;
else
    trials = detrend(trials);
end

averaged = mean(trials, 2, 'omitnan');
end

% -------------------------------------------------------------------------
function [fraction, area] = localLongestRun(mask, values, level, sign)
%LOCALLONGESTRUN Longest contiguous suprathreshold run and its excess area.

fraction = 0;
area = zeros(0, 1);
if ~any(mask)
    return
end

edges = diff([false; mask(:); false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
[runLength, which] = max(stops - starts + 1);
if isempty(which)
    return
end

fraction = runLength / numel(mask);
selected = starts(which):stops(which);
if sign > 0
    area = values(selected) - level;
else
    area = values(selected) + level;
end
end
