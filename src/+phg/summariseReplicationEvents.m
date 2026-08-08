function summary = summariseReplicationEvents(cfg, options)
%SUMMARISEREPLICATIONEVENTS Task-event density in the replication cohort.
%   The manuscript states how closely spaced the replication cohort's task
%   events are, because that spacing is the reason the inter-trial restriction
%   cannot isolate spontaneous activity there. That figure was originally typed
%   from a single run; this exports it from all of them.
%
%   Reads only the event onsets from each staged run, so it costs nothing
%   compared with the measurement.

arguments
    cfg (1,1) struct
    options.stageDir (1,1) string = cfg.replication.stageDir
    options.writeTables (1,1) logical = true
end

listing = dir(fullfile(options.stageDir, 'sub-*'));
listing = listing([listing.isdir]);

tag = strings(0, 1);
subject = strings(0, 1);
task = strings(0, 1);
nEvents = [];
medianInterval = [];
spanSeconds = [];

for k = 1:numel(listing)
    runDir = fullfile(options.stageDir, listing(k).name);
    h5File = fullfile(runDir, 'run.h5');
    metaFile = fullfile(runDir, 'meta.json');
    if ~isfile(h5File) || ~isfile(metaFile)
        continue
    end
    meta = jsondecode(fileread(metaFile));
    onsets = sort(double(h5read(h5File, '/event_onsets')));
    if numel(onsets) < 2
        continue
    end
    intervals = diff(onsets);
    % Simultaneous triggers share a timestamp and contribute zero-length
    % intervals; they are not separate opportunities for a quiet epoch, so the
    % spacing that matters is between distinct event times.
    intervals = intervals(intervals > 0);

    tag(end + 1, 1) = string(meta.tag); %#ok<AGROW>
    subject(end + 1, 1) = string(meta.subject); %#ok<AGROW>
    task(end + 1, 1) = string(meta.task); %#ok<AGROW>
    nEvents(end + 1, 1) = numel(onsets); %#ok<AGROW>
    medianInterval(end + 1, 1) = median(intervals); %#ok<AGROW>
    spanSeconds(end + 1, 1) = onsets(end) - onsets(1); %#ok<AGROW>
end

summary = table(tag, subject, task, nEvents, medianInterval, spanSeconds, ...
    'VariableNames', {'run_tag', 'subject', 'task', 'n_events', ...
    'median_inter_event_interval_seconds', 'event_span_seconds'});

if options.writeTables
    phg.writeTableAtomic(summary, ...
        fullfile(cfg.tableDir, 'replication_event_density.csv'));
end

fprintf(['[PHG] Replication event density: %d runs, median inter-event ', ...
    'interval %.2f s (range %.2f-%.2f)\n'], height(summary), ...
    median(summary.median_inter_event_interval_seconds), ...
    min(summary.median_inter_event_interval_seconds), ...
    max(summary.median_inter_event_interval_seconds));
end
