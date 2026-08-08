function outputs = runReplicationGlobalityCheck(cfg, options)
%RUNREPLICATIONGLOBALITYCHECK Is the peri-peak response contact-specific?
%
% The replication cohort produces a far higher excursion rate than the
% discovery cohort (three quarters of contacts against under a third), and in
% two participants essentially every contact is both suprathreshold and
% dilation-linked. That pattern has an obvious alternative explanation: a
% single task-evoked pupil dilation shared by the whole implant, picked up
% wherever a high-gamma peak happens to fall, rather than coupling that differs
% from contact to contact.
%
% This turns that suspicion into a measurement. Within each participant, the
% peri-peak curves of individual contacts are compared with one another. If the
% response is contact-specific, curves should differ; if it is one global
% dilation, they should be near-copies of each other and a single component
% should account for almost all of their variance.
%
% Reported for both cohorts is not possible -- the discovery archive stores
% only summary statistics per contact, not the curves -- so this is a statement
% about the replication cohort alone, and it is used to explain why the cohort
% cannot test H1 rather than as evidence about the discovery result.

arguments
    cfg (1,1) struct
    options.stageDir (1,1) string = cfg.replication.stageDir
    options.pupilVariant (1,1) string = "gaze_regressed"
    options.contactsPerSubject (1,1) double = 16
    options.writeTables (1,1) logical = true
    % Subjects can be run one process at a time. The measurement is memory
    % hungry enough that a single long process was being killed part way
    % through, losing the subjects it had already finished.
    options.subjects string = strings(0, 1)
    options.append (1,1) logical = false
end

% This diagnostic compares the *shape* of peri-peak curves, so it does not need
% the permutation count the inferential measurement uses. Twenty permutations
% give a stable mean curve at a fifth of the memory and time, which matters
% because the full setting was enough to exhaust this machine.
cfg.replication.numPermutations = 20;

runDirs = localListRuns(options.stageDir);
rng(cfg.randomSeed, 'twister');

subjects = strings(numel(runDirs), 1);
for r = 1:numel(runDirs)
    parts = split(string(runDirs(r)), filesep);
    subjects(r) = extractBefore(parts(end), "_task");
end
uniqueSubjects = unique(subjects);
if ~isempty(options.subjects)
    uniqueSubjects = intersect(uniqueSubjects, options.subjects, 'stable');
end

outputFile = fullfile(cfg.tableDir, 'replication_response_globality.csv');
rows = [];
for s = 1:numel(uniqueSubjects)
    subject = uniqueSubjects(s);
    theseRuns = runDirs(subjects == subject);

    [curves, labels, regions] = localSubjectCurves(theseRuns, cfg, ...
        options.pupilVariant, options.contactsPerSubject);
    if size(curves, 2) < 3
        continue
    end

    % Pairwise similarity between contacts' peri-peak curves.
    correlations = corr(curves, 'Rows', 'pairwise');
    offDiagonal = correlations(triu(true(size(correlations)), 1));
    offDiagonal = offDiagonal(isfinite(offDiagonal));

    % Variance explained by the leading component: one dominant component
    % means one shared response.
    centred = curves - mean(curves, 2, 'omitnan');
    centred(~isfinite(centred)) = 0;
    singularValues = svd(centred, 'econ');
    explained = singularValues .^ 2 / sum(singularValues .^ 2);

    isHippocampal = lower(regions) == "hippocampus";
    rows = [rows; {subject, size(curves, 2), sum(isHippocampal), ...
        median(offDiagonal), mean(offDiagonal > 0.9), explained(1), ...
        sum(explained(1:min(3, numel(explained))))}]; %#ok<AGROW>

    fprintf(['[PHG] %s: %d contacts, median pairwise r = %.3f, ' ...
        '%.0f%% of pairs r > 0.9, PC1 explains %.1f%%\n'], subject, ...
        size(curves, 2), median(offDiagonal), ...
        100 * mean(offDiagonal > 0.9), 100 * explained(1));
end

if isempty(rows)
    outputs = struct('summary', table());
    return
end

summary = cell2table(rows, 'VariableNames', {'subject_id', 'n_contacts', ...
    'n_hippocampal', 'median_pairwise_r', 'fraction_pairs_r_above_0p9', ...
    'pc1_variance_explained', 'pc1to3_variance_explained'});

if options.writeTables
    if options.append && isfile(outputFile)
        previous = readtable(outputFile, 'TextType', 'string');
        previous(ismember(previous.subject_id, summary.subject_id), :) = [];
        summary = sortrows([previous; summary], 'subject_id');
    end
    phg.writeTableAtomic(summary, outputFile);
end
outputs = struct('summary', summary);
end

% -------------------------------------------------------------------------
function [curves, labels, regions] = localSubjectCurves(runDirs, cfg, ...
    pupilVariant, maxContacts)
%LOCALSUBJECTCURVES Peri-peak curves for a sample of one subject's contacts.

curves = [];
labels = strings(0, 1);
regions = strings(0, 1);

% Peaks and pupil come from the first run only. One run is enough to ask
% whether contacts within a session share a response, and it keeps this
% diagnostic cheap next to the full measurement.
run = phg.loadReplicationRun(runDirs(1), 'includeRipple', false);
if ~run.anchorValid
    return
end

if pupilVariant == "gaze_regressed"
    pupil = run.pupilGazeRegressed;
else
    pupil = run.pupil;
end

frame = find(run.inAnalysisFrame);
% Always keep every hippocampal contact, then fill up to maxContacts with a
% deterministic spread across the rest of the implant.
hippocampal = frame(lower(run.channelRegion(frame)) == "hippocampus");
others = setdiff(frame, hippocampal, 'stable');
room = max(0, maxContacts - numel(hippocampal));
if numel(others) > room
    others = others(round(linspace(1, numel(others), room)));
end
selected = [hippocampal(:); others(:)];

neighbours = localShaftNeighbours(run);
nyquist = run.ieegFs / 2;

for k = 1:numel(selected)
    c = selected(k);
    if isempty(neighbours{c})
        continue
    end
    signalUv = double(run.ieeg(:, c)) - ...
        mean(double(run.ieeg(:, neighbours{c})), 2);
    signalUv(~isfinite(signalUv)) = 0;
    for f = cfg.replication.notchHz
        if f >= nyquist
            continue
        end
        [b, a] = iirnotch(f / nyquist, cfg.replication.notchBandwidth);
        signalUv = filtfilt(b, a, signalUv);
    end
    [b, a] = butter(cfg.replication.filterOrder, ...
        cfg.replication.highGammaHz / nyquist, 'bandpass');
    power = abs(hilbert(filtfilt(b, a, signalUv))) .^ 2;

    inEpoch = power(run.coreMask);
    if numel(inEpoch) < run.ieegFs
        continue
    end
    threshold = median(inEpoch) + iqr(inEpoch) * cfg.replication.thresholdIqr;
    above = power > threshold & run.coreMask;
    edges = diff([false; above; false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    times = nan(numel(starts), 1);
    for m = 1:numel(starts)
        [~, offset] = max(power(starts(m):stops(m)));
        times(m) = run.ieegTime(starts(m) + offset - 1);
    end

    [~, diagnostics] = phg.measurePeriPeakResponse(pupil, run.pupilTime, ...
        times, cfg);
    if isempty(diagnostics.pd)
        continue
    end
    curves = [curves, diagnostics.pd]; %#ok<AGROW>
    labels(end + 1, 1) = run.channelName(c); %#ok<AGROW>
    regions(end + 1, 1) = run.channelRegion(c); %#ok<AGROW>
    clear signalUv power filtered
end
end

% -------------------------------------------------------------------------
function neighbours = localShaftNeighbours(run)
nChannels = numel(run.channelName);
neighbours = cell(nChannels, 1);
for c = 1:nChannels
    sameShaft = find(run.channelShaft == run.channelShaft(c));
    neighbours{c} = sameShaft( ...
        abs(run.channelIndex(sameShaft) - run.channelIndex(c)) == 1);
end
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
