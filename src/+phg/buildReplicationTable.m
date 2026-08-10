function [derived, runSummary] = buildReplicationTable(cfg, options)
%BUILDREPLICATIONTABLE Measure the EBRAINS cohort with the discovery pipeline.
%
% Produces a table in the schema `runDerivedAnalyses` consumes, so the hurdle
% model runs on the replication cohort unchanged. The signal chain reproduces
% the discovery cohort step for step, verified against the stored analysis
% object (obj.Fs = 1000, obj.BandPass = [70 170], obj.RRType = 'Laplace',
% obj.ThrX = 5, obj.TimeRng = [-20000 20000]):
%
%   1. shaft Laplacian: contact minus the mean of its immediate neighbours
%   2. notch at 60, 120, 180 and 240 Hz
%   3. bandpass butter(1, [70 170]/500) applied with filtfilt
%   4. high-gamma power as abs(hilbert(.))^2
%   5. peak times as the maximum within each run above median + 5 x IQR
%   6. peri-peak pupil measurement, see measurePeriPeakResponse
%
% Contacts are measured once per run and then pooled across runs within a
% participant by concatenating peak times on a common pupil timeline, because a
% contact is one electrode regardless of how many runs it was recorded in.

arguments
    cfg (1,1) struct
    options.stageDir (1,1) string = cfg.replication.stageDir
    options.pupilVariant (1,1) string = "gaze_regressed"
    options.epochRule (1,1) string = "A"
    options.writeTables (1,1) logical = true
end

runDirs = localListRuns(options.stageDir);
if isempty(runDirs)
    % Distinguish "the drive is not mounted" from "the data is missing".
    % The staged cohort lives on external storage, so the common failure is a
    % disconnected volume rather than a lost dataset, and the two need
    % different responses.
    parent = fileparts(options.stageDir);
    if ~isfolder(parent)
        error('phg:StagingVolumeUnavailable', ...
            ['Staging volume is not mounted: %s does not exist. The ' ...
             'replication cohort lives on external storage; reconnect the ' ...
             'drive, or re-stage with tools/ingest_ebrains.py. Replication ' ...
             'results already computed are in results/tables and need no ' ...
             'access to this volume.'], parent);
    end
    error('phg:NoStagedRuns', ...
        'No staged replication runs under %s', options.stageDir);
end

fprintf('[PHG] Replication measurement: %d runs, pupil=%s, epoch rule=%s\n', ...
    numel(runDirs), options.pupilVariant, options.epochRule);

% One pupil trace is stored per run, not per contact. Every contact recorded in
% a run shares that run's pupil signal, so keeping a copy inside each contact's
% entry multiplied the pupil data by the contact count -- around a hundredfold,
% which is enough to exhaust memory on the full cohort.
accumulator = containers.Map('KeyType', 'char', 'ValueType', 'any');
runPupil = cell(numel(runDirs), 1);
runPupilTime = cell(numel(runDirs), 1);
runEvents = cell(numel(runDirs), 1);
runSubject = strings(numel(runDirs), 1);
runRows = [];

for r = 1:numel(runDirs)
    run = phg.loadReplicationRun(runDirs(r), 'includeRipple', false);
    if ~run.anchorValid
        warning('phg:ReplicationAnchor', ...
            'Run %s failed time-anchor validation and is skipped.', run.tag);
        continue
    end

    switch options.epochRule
        case "A", mask = run.coreMask;
        case "B", mask = run.maskRuleB;
        case "C", mask = run.maskRuleC;
        otherwise
            error('phg:UnknownEpochRule', ...
                'Epoch rule must be A, B or C; got %s', options.epochRule);
    end

    if options.pupilVariant == "gaze_regressed"
        pupil = run.pupilGazeRegressed;
    else
        pupil = run.pupil;
    end

    [peakTimes, envelopeStats] = localDetectPeaks(run, mask, cfg);

    runPupil{r} = pupil;
    runPupilTime{r} = run.pupilTime;
    runEvents{r} = run.eventOnsets;
    runSubject(r) = run.subject;

    for c = 1:numel(run.channelName)
        if ~run.inAnalysisFrame(c)
            continue
        end
        key = char(run.subject + "|" + run.channelName(c));
        entry = struct('subject', run.subject, ...
            'channel', run.channelName(c), 'region', run.channelRegion(c), ...
            'shaft', run.channelShaft(c), 'xyz', run.channelXyz(c, :), ...
            'peakTimes', [], 'runIndices', [], 'runTag', "");
        if accumulator.isKey(key)
            entry = accumulator(key);
        end
        % Peak times stay in their own run's time base and are remapped onto
        % the stitched axis at measurement time. The contact records which runs
        % it appears in rather than carrying their signals around.
        entry.peakTimes = [entry.peakTimes; {peakTimes{c}}];
        entry.runIndices = [entry.runIndices; r];
        entry.runTag = entry.runTag + run.tag + ";";
        accumulator(key) = entry;
    end

    % How much of the pupil signal is gaze? The discovery cohort lists
    % uncontrolled gaze as a limitation and could not answer this. Reported as
    % variance removed by the regression, whatever the answer turns out to be.
    both = isfinite(run.pupil) & isfinite(run.pupilGazeRegressed);
    gazeVarianceRemoved = NaN;
    if sum(both) > 100
        rawVariance = var(run.pupil(both));
        if rawVariance > 0
            gazeVarianceRemoved = max(0, 1 - ...
                var(run.pupilGazeRegressed(both)) / rawVariance);
        end
    end

    runRows = [runRows; {run.tag, run.subject, run.task, run.runIndex, ...
        sum(mask) / run.ieegFs, numel(run.channelName), ...
        sum(run.inAnalysisFrame), median(envelopeStats.nPeaks), ...
        mean(isnan(pupil)), gazeVarianceRemoved}]; %#ok<AGROW>

    % The wideband matrix for a run is hundreds of megabytes and nothing below
    % needs it; release it before the next run is loaded.
    clear run pupil mask
end

runSummary = cell2table(runRows, 'VariableNames', {'run_tag', 'subject', ...
    'task', 'run_index', 'analysed_seconds', 'staged_channels', ...
    'analysis_contacts', 'median_peaks_per_contact', ...
    'pupil_missing_fraction', 'gaze_variance_removed'});

derived = localMeasureContacts(accumulator, runPupil, runPupilTime, ...
    runEvents, runSubject, cfg);

if options.writeTables
    suffix = "_" + options.pupilVariant + "_rule" + options.epochRule;
    phg.writeTableAtomic(derived, fullfile(cfg.tableDir, ...
        "replication_contact_measures" + suffix + ".csv"));
    phg.writeTableAtomic(runSummary, fullfile(cfg.tableDir, ...
        'replication_run_summary.csv'));
end
end

% -------------------------------------------------------------------------
function [peakTimes, stats] = localDetectPeaks(run, mask, cfg)
%LOCALDETECTPEAKS Discovery signal chain, then median + 5 x IQR peak times.

fs = run.ieegFs;
nChannels = numel(run.channelName);
peakTimes = cell(nChannels, 1);
peakCount = zeros(nChannels, 1);

nyquist = fs / 2;
notchFrequencies = cfg.replication.notchHz;
band = cfg.replication.highGammaHz;

% The Laplacian is formed one channel at a time. Materialising it for the whole
% run allocates a double-precision copy of a matrix that is already hundreds of
% megabytes in single precision, which is enough to end the job.
neighbours = localShaftNeighbours(run);

for c = 1:nChannels
    if isempty(neighbours{c})
        peakTimes{c} = [];
        continue
    end
    signalUv = double(run.ieeg(:, c)) - ...
        mean(double(run.ieeg(:, neighbours{c})), 2);
    if all(~isfinite(signalUv))
        peakTimes{c} = [];
        continue
    end
    signalUv(~isfinite(signalUv)) = 0;

    for f = notchFrequencies
        if f >= nyquist
            continue
        end
        [b, a] = iirnotch(f / nyquist, cfg.replication.notchBandwidth);
        signalUv = filtfilt(b, a, signalUv);
    end

    [b, a] = butter(cfg.replication.filterOrder, band / nyquist, 'bandpass');
    filtered = filtfilt(b, a, signalUv);
    power = abs(hilbert(filtered)) .^ 2;

    % The threshold and the peaks both come from eligible epochs only: padding
    % exists to give the filters context, not to contribute observations.
    inEpoch = power(mask);
    if numel(inEpoch) < fs
        peakTimes{c} = [];
        continue
    end
    threshold = median(inEpoch) + iqr(inEpoch) * cfg.replication.thresholdIqr;

    above = power > threshold & mask;
    edges = diff([false; above; false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;

    times = nan(numel(starts), 1);
    for k = 1:numel(starts)
        [~, offset] = max(power(starts(k):stops(k)));
        times(k) = run.ieegTime(starts(k) + offset - 1);
    end
    peakTimes{c} = times;
    peakCount(c) = numel(times);
end

stats = struct('nPeaks', peakCount);
end

% -------------------------------------------------------------------------
function neighbours = localShaftNeighbours(run)
%LOCALSHAFTNEIGHBOURS Immediate neighbours of each contact on its own shaft.
%   Mirrors PupilHG.m:323-345: two neighbours where both exist, the single
%   available neighbour at a shaft end, and none where neither exists -- in
%   which case the contact has no Laplacian and is dropped, exactly as the
%   discovery code leaves it NaN.

nChannels = numel(run.channelName);
neighbours = cell(nChannels, 1);
for c = 1:nChannels
    sameShaft = find(run.channelShaft == run.channelShaft(c));
    here = run.channelIndex(c);
    neighbours{c} = sameShaft( ...
        abs(run.channelIndex(sameShaft) - here) == 1);
end
end

% -------------------------------------------------------------------------
function derived = localMeasureContacts(accumulator, runPupil, runPupilTime, ...
    runEvents, runSubject, cfg)
%LOCALMEASURECONTACTS Run the peri-peak measurement for every contact.
%   Contacts belonging to one participant share that participant's runs, so the
%   stitched pupil timeline is built once per participant and reused, rather
%   than rebuilt for each of their ~80 contacts.

keys = accumulator.keys;
n = numel(keys);
PtID = strings(n, 1);
Label = strings(n, 1);
NMM = strings(n, 1);
Shaft = strings(n, 1);
XYZMNI = nan(n, 3);
RespSig = nan(n, 1);
RespAreaNet = nan(n, 1);
RespAreaAbs = nan(n, 1);
NPeaks = nan(n, 1);
NPeaksUsed = nan(n, 1);
MeanTrialCount = nan(n, 1);
WindowMissingFraction = nan(n, 1);
EventOverlap = nan(n, 1);
FitHeight = nan(n, 1);
FitR2 = nan(n, 1);

% Order contacts by participant so each stitched timeline is built once.
subjectOf = strings(n, 1);
for k = 1:n
    subjectOf(k) = accumulator(keys{k}).subject;
end
[~, order] = sort(subjectOf);
cachedSubject = "";
cachedPupil = [];
cachedTime = [];
cachedOffsets = [];
cachedRuns = [];

for index = 1:n
    k = order(index);
    entry = accumulator(keys{k});
    PtID(k) = entry.subject;
    Label(k) = entry.channel;
    NMM(k) = entry.region;
    Shaft(k) = entry.shaft;
    XYZMNI(k, :) = entry.xyz;

    if entry.subject ~= cachedSubject
        cachedRuns = find(runSubject == entry.subject);
        [cachedPupil, cachedTime, cachedOffsets] = localStitchRuns( ...
            runPupil(cachedRuns), runPupilTime(cachedRuns), cfg);
        cachedSubject = entry.subject;
    end

    peaks = localRemapPeaks(entry.peakTimes, entry.runIndices, cachedRuns, ...
        cachedOffsets, runPupilTime, cfg);

    [response, ~] = phg.measurePeriPeakResponse(cachedPupil, cachedTime, ...
        peaks, cfg);
    EventOverlap(k) = localEventOverlap(entry.peakTimes, ...
        runEvents(entry.runIndices), cfg);
    FitHeight(k) = response.FitHeight;
    FitR2(k) = response.FitR2;
    RespSig(k) = response.RespSig;
    RespAreaNet(k) = response.RespAreaNet;
    RespAreaAbs(k) = response.RespAreaAbs;
    NPeaks(k) = response.NPeaks;
    NPeaksUsed(k) = response.NPeaksUsed;
    MeanTrialCount(k) = response.MeanTrialCount;
    WindowMissingFraction(k) = response.WindowMissingFraction;

    if mod(index, 50) == 0
        fprintf('[PHG]   measured %d/%d contacts\n', index, n);
    end
end

derived = table(PtID, Label, NMM, Shaft, XYZMNI, RespSig, RespAreaNet, ...
    RespAreaAbs, NPeaks, NPeaksUsed, MeanTrialCount, WindowMissingFraction, ...
    EventOverlap, FitHeight, FitR2);
derived.Chan = (1:height(derived))';
derived.FitShift = nan(height(derived), 1);
derived = sortrows(derived, {'PtID', 'Label'});
end

% -------------------------------------------------------------------------
function overlap = localEventOverlap(peakCells, eventCells, cfg)
%LOCALEVENTOVERLAP Share of the peri-peak window sitting near a task event.
%   Required by docs/replication_amendment_01.md section 3: because the +/-2 s
%   restriction applies to peak times rather than to whole windows, the peri-peak
%   average can overlap task-evoked pupil responses. The amount of overlap is
%   measured so that a difference between contact classes cannot go unnoticed,
%   and enters the adjusted model if it exceeds the prespecified 5-point gap.
%
%   Computed within each run against that run's own events, then pooled across
%   the contact's peaks. Doing it on the stitched timeline would compare peaks
%   against events from a different recording.

guard = 2.0;
halfWindow = cfg.replication.responseWindowSeconds;
sampleLags = linspace(-halfWindow, halfWindow, 101);

allFractions = [];
for run = 1:numel(peakCells)
    peakTimes = peakCells{run};
    eventOnsets = eventCells{run};
    if isempty(peakTimes) || isempty(eventOnsets)
        continue
    end

    % Sub-sample when a run contributes very many peaks; this is a mean over
    % peaks and converges long before the full set is needed.
    if numel(peakTimes) > 300
        peakTimes = peakTimes(round(linspace(1, numel(peakTimes), 300)));
    end

    % Event onsets are not unique -- simultaneous triggers share a timestamp --
    % so nearest-neighbour lookup is by binary search rather than interp1,
    % which rejects duplicated sample points outright.
    ordered = unique(sort(eventOnsets));
    if numel(ordered) < 2
        allFractions = [allFractions; ...
            double(abs(peakTimes(:) - ordered(1)) <= guard)]; %#ok<AGROW>
        continue
    end

    fractions = nan(numel(peakTimes), 1);
    for k = 1:numel(peakTimes)
        query = peakTimes(k) + sampleLags;
        position = min(max(discretize(query, [-inf; ordered; inf]) - 1, 1), ...
            numel(ordered) - 1);
        distance = min(abs(query - ordered(position)'), ...
            abs(query - ordered(position + 1)'));
        fractions(k) = mean(distance <= guard);
    end
    allFractions = [allFractions; fractions]; %#ok<AGROW>
end

if isempty(allFractions)
    overlap = NaN;
else
    overlap = mean(allFractions, 'omitnan');
end
end

% -------------------------------------------------------------------------
function [pupil, pupilTime, startSample] = localStitchRuns(pupilCells, ...
    timeCells, cfg)
%LOCALSTITCHRUNS Lay a participant's runs end to end on one uniform timeline.
%   A contact is one electrode however many runs it was recorded in, and the
%   discovery measurement averages over all of its peaks at once, so the runs
%   have to share a timeline. They are laid end to end separated by a guard
%   band of NaN wider than the full peri-peak window.
%
%   The guard is what makes this safe: any window that would straddle a join
%   reaches into NaN and is rejected by the missing-data check in
%   measurePeriPeakResponse rather than splicing the end of one run onto the
%   start of another.
%
%   STARTSAMPLE gives the zero-based sample index at which each run begins on
%   the stitched axis, which is what localRemapPeaks needs to place peaks.

fs = cfg.replication.pupilFs;
guardSeconds = 2 * abs(cfg.replication.timeRangeMs(1)) / 1000 + 5;
guardSamples = round(guardSeconds * fs);

nRuns = numel(pupilCells);
startSample = zeros(nRuns, 1);
segments = cell(nRuns, 1);
cursor = 0;
for k = 1:nRuns
    values = pupilCells{k}(:);
    startSample(k) = cursor;
    cursor = cursor + numel(values);
    if k < nRuns
        values = [values; nan(guardSamples, 1)]; %#ok<AGROW>
        cursor = cursor + guardSamples;
    end
    segments{k} = values;
end

pupil = vertcat(segments{:});
pupilTime = (0:numel(pupil) - 1)' / fs;
if nRuns == 1
    % Single run: keep its native time base so peak times need no shifting.
    pupilTime = timeCells{1}(:);
end
end

% -------------------------------------------------------------------------
function peaks = localRemapPeaks(peakCells, runIndices, subjectRuns, ...
    startSample, runPupilTime, cfg)
%LOCALREMAPPEAKS Place a contact's peak times on the stitched timeline.
%   Peak times are recorded against each run's own pupil clock. On the stitched
%   axis a peak sits at its offset into its own run, plus the number of samples
%   laid down before that run began.

fs = cfg.replication.pupilFs;
peaks = [];
for k = 1:numel(runIndices)
    times = peakCells{k}(:);
    if isempty(times)
        continue
    end
    slot = find(subjectRuns == runIndices(k), 1);
    if isempty(slot)
        continue
    end
    if isscalar(subjectRuns)
        peaks = [peaks; times]; %#ok<AGROW>
        continue
    end
    runStart = runPupilTime{runIndices(k)}(1);
    peaks = [peaks; (times - runStart) + startSample(slot) / fs]; %#ok<AGROW>
end
peaks = sort(peaks);
end

% -------------------------------------------------------------------------
function runDirs = localListRuns(stageDir)
%LOCALLISTRUNS Staged run directories, in deterministic order.

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
