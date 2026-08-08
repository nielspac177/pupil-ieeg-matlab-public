function outputs = runRipplePupilTest(cfg, options)
%RUNRIPPLEPUPILTEST H2: are hippocampal ripples enriched during constriction?
%
% The mechanism proposed in the discovery Discussion, tested directly for the
% first time. The discovery archive held trial-averaged spectra, so the only
% available proxy was whether hippocampal coupling was spectrally concentrated
% near the ripple band; that came out null, and averaging over hundreds of
% high-gamma peaks blurs exactly the band-limited structure such a test looks
% for. This cohort has wideband signal, so ripples are detectable as discrete
% events and the question can be asked as it was originally posed.
%
% Test, prespecified in docs/replication_plan_ebrains.md section 2:
%   For every detected hippocampal ripple, take the pupil derivative in a
%   window around it. Compare the fraction occurring during constriction
%   against a null built by circularly shifting ripple times within the run.
%
% Secondary and labelled as such regardless of outcome. A null here means the
% ripple paragraph in the Discussion should be cut, not softened.
%
% Referencing, fixed before any ripple was detected. Primary detection runs on
% the shaft Laplacian at 1 kHz: it matches the referencing used for the
% high-gamma analysis, and a common reference is a poor basis for ripple
% detection because volume-conducted and myogenic transients appear on every
% contact at once. 1 kHz resolves the 80-120 Hz band with room to spare. The
% 4 kHz common-reference signal staged for hippocampal contacts is run as a
% sensitivity analysis and reported alongside.

arguments
    cfg (1,1) struct
    options.stageDir (1,1) string = cfg.replication.stageDir
    options.pupilVariant (1,1) string = "gaze_regressed"
    options.signalSource (1,1) string {mustBeMember(options.signalSource, ...
        ["laplacian_1khz", "raw_4khz"])} = "laplacian_1khz"
    options.writeTables (1,1) logical = true
end

label = options.signalSource;
runDirs = localListRuns(options.stageDir);
rng(cfg.randomSeed, 'twister');

eventRows = [];
runRows = [];
perRun = struct('rippleTimes', {}, 'pupilTime', {}, 'derivative', {}, ...
    'span', {}, 'origin', {});
% phg.detectRipples had never been run on real data before this analysis, so
% its behaviour is reported rather than assumed: every candidate carries an
% explicit rejection reason, and the tally is exported.
rejectionTally = containers.Map('KeyType', 'char', 'ValueType', 'double');

for r = 1:numel(runDirs)
    run = phg.loadReplicationRun(runDirs(r));
    if isempty(run.rippleChannels) || ~run.anchorValid
        continue
    end

    [signalUv, fs, timeAxis, channelNames] = ...
        localRippleSignal(run, options.signalSource);
    if isempty(signalUv)
        continue
    end

    if options.pupilVariant == "gaze_regressed"
        pupil = run.pupilGazeRegressed;
    else
        pupil = run.pupil;
    end

    % Pupil derivative, smoothed over 500 ms before differencing so that the
    % sign reports the direction of the slow pupil movement rather than sample
    % noise. Constriction is a negative derivative.
    smoothed = movmean(pupil, round(0.5 * run.pupilFs), 'omitnan');
    derivative = [NaN; diff(smoothed)] * run.pupilFs;

    % Detection is restricted to eligible epochs, matched to the high-gamma
    % analysis, and is never restricted by any pupil criterion.
    epochMask = interp1(run.ieegTime, double(run.coreMask), timeAxis, ...
        'nearest', 0) > 0.5;

    events = phg.detectRipples(signalUv, fs, ...
        'bandHz', cfg.ripple.bandHz, ...
        'filterOrder', cfg.ripple.filterOrder, ...
        'rmsWindowSeconds', cfg.ripple.rmsWindowSeconds, ...
        'minimumZ', cfg.ripple.minimumZ, ...
        'maximumZ', cfg.ripple.maximumZ, ...
        'minimumDurationSeconds', cfg.ripple.minimumDurationSeconds, ...
        'maximumRawAmplitudeUv', cfg.ripple.maximumRawAmplitudeUv, ...
        'spectralWindowSeconds', cfg.ripple.spectralWindowSeconds);

    reasons = string(events.rejection_reason);
    for k = 1:numel(reasons)
        key = char(reasons(k));
        if rejectionTally.isKey(key)
            rejectionTally(key) = rejectionTally(key) + 1;
        else
            rejectionTally(key) = 1;
        end
    end

    accepted = events(events.accepted, :);
    if isempty(accepted)
        runRows = [runRows; {run.tag, run.subject, run.task, ...
            numel(channelNames), height(events), 0, NaN, NaN}]; %#ok<AGROW>
        continue
    end

    centre = accepted.onset_sample + ...
        round((accepted.offset_sample - accepted.onset_sample) / 2);
    centre = min(max(centre, 1), numel(timeAxis));
    inEpoch = epochMask(centre);
    accepted = accepted(inEpoch, :);
    if isempty(accepted)
        runRows = [runRows; {run.tag, run.subject, run.task, ...
            numel(channelNames), height(events), 0, NaN, NaN}]; %#ok<AGROW>
        continue
    end
    rippleTimes = timeAxis(centre(inEpoch));

    derivativeAtRipple = interp1(run.pupilTime, derivative, rippleTimes, ...
        'linear', NaN);
    valid = isfinite(derivativeAtRipple);
    observedFraction = mean(derivativeAtRipple(valid) < 0);

    eventRows = [eventRows; table( ...
        repmat(run.subject, height(accepted), 1), ...
        repmat(run.tag, height(accepted), 1), ...
        channelNames(accepted.channel), rippleTimes, ...
        accepted.duration_seconds, accepted.peak_z, ...
        accepted.peak_frequency_hz, derivativeAtRipple, ...
        'VariableNames', {'subject', 'run_tag', 'channel', 'ripple_time', ...
        'duration_seconds', 'peak_z', 'peak_frequency_hz', ...
        'pupil_derivative'})]; %#ok<AGROW>

    % Kept so that the pooled null can shift ripple times in real time rather
    % than permuting the derivative values that were sampled at them.
    perRun(end + 1) = struct('rippleTimes', rippleTimes, ...
        'pupilTime', run.pupilTime, 'derivative', derivative, ...
        'span', max(timeAxis) - min(timeAxis), ...
        'origin', min(timeAxis)); %#ok<AGROW>

    runRows = [runRows; {run.tag, run.subject, run.task, ...
        numel(channelNames), height(events), height(accepted), ...
        observedFraction, NaN}]; %#ok<AGROW>
end

if isempty(eventRows)
    warning('phg:NoRipples', ...
        'No accepted ripples in eligible epochs; H2 is not estimable.');
    outputs = struct('events', table(), 'runSummary', table(), ...
        'test', table(), 'nullDistribution', [], 'label', label);
    return
end

runSummary = cell2table(runRows, 'VariableNames', {'run_tag', 'subject', ...
    'task', 'n_hippocampal_channels', 'n_candidate_events', ...
    'n_accepted_in_epoch', 'constriction_fraction', 'null_mean_fraction'});

valid = isfinite(eventRows.pupil_derivative);
observed = mean(eventRows.pupil_derivative(valid) < 0);

% Pooled circular-shift null. Every ripple time in a run is displaced by one
% common random offset and wrapped within that run, which preserves the ripple
% rate and the clustering of ripples in time while destroying their alignment
% to the pupil. Counts are pooled across runs so the surrogate statistic is
% built the same way as the observed one.
nShift = cfg.replication.numCircularShifts;
pooledNull = nan(nShift, 1);
for s = 1:nShift
    hits = 0;
    total = 0;
    for k = 1:numel(perRun)
        entry = perRun(k);
        if entry.span <= 0
            continue
        end
        shifted = entry.origin + ...
            mod(entry.rippleTimes - entry.origin + rand * entry.span, ...
            entry.span);
        sampled = interp1(entry.pupilTime, entry.derivative, shifted, ...
            'linear', NaN);
        usable = isfinite(sampled);
        hits = hits + sum(sampled(usable) < 0);
        total = total + sum(usable);
    end
    if total > 0
        pooledNull(s) = hits / total;
    end
end

% One-sided: the hypothesis is enrichment during constriction.
pValue = (1 + sum(pooledNull >= observed)) / (1 + sum(isfinite(pooledNull)));

test = table(string(label), "hippocampal_ripples_during_constriction", ...
    sum(valid), observed, mean(pooledNull, 'omitnan'), ...
    std(pooledNull, 'omitnan'), pValue, nShift, ...
    'VariableNames', {'signal_source', 'test', 'n_ripples', ...
    'observed_constriction_fraction', 'null_mean_fraction', ...
    'null_sd_fraction', 'p_value_one_sided', 'n_circular_shifts'});

if options.writeTables
    suffix = "_" + label;
    phg.writeTableAtomic(runSummary, fullfile(cfg.tableDir, ...
        "replication_ripple_run_summary" + suffix + ".csv"));
    phg.writeTableAtomic(test, fullfile(cfg.tableDir, ...
        "replication_ripple_pupil_test" + suffix + ".csv"));
    % The surrogate distribution is exported in full: the figure plots it, and
    % a reader can check the p-value against the observed statistic.
    nullTable = table((1:numel(pooledNull))', pooledNull, ...
        'VariableNames', {'surrogate_index', 'constriction_fraction'});
    phg.writeTableAtomic(nullTable, fullfile(cfg.tableDir, ...
        "replication_ripple_null_distribution" + suffix + ".csv"));

    reasonKeys = string(rejectionTally.keys)';
    reasonCounts = cell2mat(rejectionTally.values)';
    [reasonCounts, order] = sort(reasonCounts, 'descend');
    detectorTable = table(reasonKeys(order), reasonCounts, ...
        reasonCounts / sum(reasonCounts), 'VariableNames', ...
        {'outcome', 'n_candidates', 'proportion_of_candidates'});
    phg.writeTableAtomic(detectorTable, fullfile(cfg.tableDir, ...
        "replication_ripple_detector_audit" + suffix + ".csv"));
end

fprintf(['[PHG] H2 (%s): %d ripples, %.1f%% during constriction against a ', ...
    'null of %.1f%%, shift-null P = %.3g\n'], label, sum(valid), ...
    100 * observed, 100 * mean(pooledNull, 'omitnan'), pValue);

outputs = struct('events', eventRows, 'runSummary', runSummary, ...
    'test', test, 'nullDistribution', pooledNull, 'label', label);
end

% -------------------------------------------------------------------------
function [signalUv, fs, timeAxis, channelNames] = localRippleSignal(run, source)
%LOCALRIPPLESIGNAL Hippocampal signal for ripple detection, in one reference.

channelNames = run.rippleChannels;
if isempty(channelNames)
    signalUv = [];
    fs = NaN;
    timeAxis = [];
    return
end

if source == "raw_4khz"
    signalUv = double(run.rippleIeeg);
    fs = run.rippleFs;
    timeAxis = run.rippleTime;
    return
end

% Shaft Laplacian at 1 kHz, formed the same way as for the high-gamma
% analysis: contact minus the mean of its immediate neighbours on the shaft.
fs = run.ieegFs;
timeAxis = run.ieegTime;
signalUv = nan(numel(timeAxis), numel(channelNames));
keep = true(numel(channelNames), 1);
for k = 1:numel(channelNames)
    c = find(run.channelName == channelNames(k), 1);
    if isempty(c)
        keep(k) = false;
        continue
    end
    sameShaft = find(run.channelShaft == run.channelShaft(c));
    neighbours = sameShaft( ...
        abs(run.channelIndex(sameShaft) - run.channelIndex(c)) == 1);
    if isempty(neighbours)
        keep(k) = false;
        continue
    end
    signalUv(:, k) = double(run.ieeg(:, c)) - ...
        mean(double(run.ieeg(:, neighbours)), 2);
end
signalUv = signalUv(:, keep);
channelNames = channelNames(keep);
end

% -------------------------------------------------------------------------
function runDirs = localListRuns(stageDir)
listing = dir(fullfile(stageDir, 'sub-*'));
listing = listing([listing.isdir]);
runDirs = strings(0, 1);
for k = 1:numel(listing)
    candidate = fullfile(stageDir, listing(k).name);
    if isfile(fullfile(candidate, 'run.h5'))
        runDirs(end + 1, 1) = string(candidate); %#ok<AGROW>
    end
end
runDirs = sort(runDirs);
end
