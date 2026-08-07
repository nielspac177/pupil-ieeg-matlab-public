function files = exportBidsRun(runData, cfg)
%EXPORTBIDSRUN Export one de-identified iEEG+pupil run to BIDS staging.
% This function writes EEGLAB SET/FDT plus BIDS sidecars. It requires that
% identifiers, dates, event messages, and anatomy have already passed the
% project's governance and de-identification gates.

required = ["subject", "session", "task", "dataUv", "samplingFrequency", ...
    "channelNames", "channelTypes", "electrodes", "events", "eyeTracking"];
missing = required(~isfield(runData, cellstr(required)));
if ~isempty(missing)
    error('phg:InvalidBidsRun', 'BIDS run lacks: %s', strjoin(missing, ', '));
end
if exist('pop_importdata', 'file') ~= 2 || exist('pop_saveset', 'file') ~= 2
    error('phg:MissingEEGLAB', ...
        'EEGLAB must be on the MATLAB path for SET/FDT export.');
end

subject = normalizeLabel(runData.subject, 'sub-');
session = normalizeLabel(runData.session, 'ses-');
task = regexprep(lower(string(runData.task)), '[^a-zA-Z0-9+]', '');
if task == ""
    error('phg:InvalidBidsLabel', 'Task label is empty after normalization.');
end

runDir = fullfile(cfg.bidsRoot, subject, session, 'ieeg');
if ~isfolder(runDir)
    mkdir(runDir);
end
stem = subject + "_" + session + "_task-" + task;

data = double(runData.dataUv);
if height(data) < width(data)
    warning('phg:DataOrientation', ...
        'Expected samples-by-channels; input has fewer rows than columns.');
end
if width(data) ~= numel(runData.channelNames)
    error('phg:ChannelCountMismatch', ...
        'dataUv columns must match channelNames.');
end

EEG = pop_importdata('dataformat', 'array', 'data', data', ...
    'srate', runData.samplingFrequency, 'setname', char(stem));
for k = 1:numel(runData.channelNames)
    EEG.chanlocs(k).labels = char(string(runData.channelNames(k)));
end
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG, 'filename', char(stem + "_ieeg.set"), ...
    'filepath', char(runDir), 'savemode', 'twofiles');

setPath = fullfile(runDir, stem + "_ieeg.set");
jsonPath = fullfile(runDir, stem + "_ieeg.json");
channelsPath = fullfile(runDir, stem + "_channels.tsv");
eventsPath = fullfile(runDir, stem + "_events.tsv");
electrodePath = fullfile(runDir, stem + "_space-MNI152Lin_electrodes.tsv");
coordinatePath = fullfile(runDir, stem + "_space-MNI152Lin_coordsystem.json");
eyePath = fullfile(runDir, stem + "_recording-eye1_physio.tsv.gz");
eyeJsonPath = fullfile(runDir, stem + "_recording-eye1_physio.json");

ieegSidecar = struct;
ieegSidecar.TaskName = char(runData.task);
ieegSidecar.SamplingFrequency = runData.samplingFrequency;
ieegSidecar.PowerLineFrequency = getOr(runData, 'powerLineFrequency', 60);
ieegSidecar.SoftwareFilters = 'n/a';
ieegSidecar.iEEGReference = char(getOr(runData, 'reference', 'unknown'));
ieegSidecar.RecordingDuration = height(data) / runData.samplingFrequency;
ieegSidecar.RecordingType = 'continuous';
ieegSidecar.Manufacturer = char(getOr(runData, 'manufacturer', 'Blackrock'));
ieegSidecar.ECOGChannelCount = sum(contains(upper(string(runData.channelTypes)), 'ECOG'));
ieegSidecar.SEEGChannelCount = sum(contains(upper(string(runData.channelTypes)), 'SEEG'));
phg.writeJson(jsonPath, ieegSidecar, true);

channels = table(string(runData.channelNames(:)), ...
    upper(string(runData.channelTypes(:))), ...
    repmat("uV", numel(runData.channelNames), 1), ...
    repmat(runData.samplingFrequency, numel(runData.channelNames), 1), ...
    repmat("good", numel(runData.channelNames), 1), ...
    strings(numel(runData.channelNames), 1), ...
    'VariableNames', {'name', 'type', 'units', 'sampling_frequency', ...
    'status', 'status_description'});
if isfield(runData, 'badChannelMask')
    channels.status(runData.badChannelMask) = "bad";
    channels.status_description(runData.badChannelMask) = "excluded by QC";
end
writetable(channels, channelsPath, 'FileType', 'text', 'Delimiter', '\t');

events = runData.events;
assertColumns(events, ["onset", "duration"]);
writetable(events, eventsPath, 'FileType', 'text', 'Delimiter', '\t');

electrodes = runData.electrodes;
assertColumns(electrodes, ["name", "x", "y", "z", "size"]);
writetable(electrodes, electrodePath, 'FileType', 'text', 'Delimiter', '\t');
coordinateSidecar = struct( ...
    'iEEGCoordinateSystem', 'MNI152Lin', ...
    'iEEGCoordinateUnits', 'mm', ...
    'iEEGCoordinateSystemDescription', ...
    'MNI coordinates exported from the laboratory localization workflow.', ...
    'iEEGCoordinateProcessingDescription', ...
    'See electrode localization provenance; projected contacts are labeled in derivatives.', ...
    'iEEGCoordinateProcessingReference', ...
    'LeGUI: Frontiers in Neuroscience 2021;15:769872');
phg.writeJson(coordinatePath, coordinateSidecar, true);

eye = runData.eyeTracking;
assertColumns(eye, ["timestamp", "x_coordinate", "y_coordinate", "pupil_size"]);
temporaryEye = fullfile(runDir, stem + "_recording-eye1_physio.tsv");
writetable(eye, temporaryEye, 'FileType', 'text', 'Delimiter', '\t');
gzip(temporaryEye, runDir);
delete(temporaryEye);
eyeSidecar = struct;
eyeSidecar.PhysioType = 'eyetrack';
eyeSidecar.RecordedEye = char(getOr(runData, 'recordedEye', 'left'));
eyeSidecar.SamplingFrequency = getOr(runData, 'eyeSamplingFrequency', ...
    1 / median(diff(eye.timestamp), 'omitnan'));
eyeSidecar.StartTime = eye.timestamp(1);
eyeSidecar.Columns = {'timestamp', 'x_coordinate', 'y_coordinate', 'pupil_size'};
eyeSidecar.x_coordinate = struct('Units', 'pixels');
eyeSidecar.y_coordinate = struct('Units', 'pixels');
eyeSidecar.pupil_size = struct('Description', ...
    'Pupil size; area versus diameter must be specified before release.', ...
    'Units', 'arbitrary');
phg.writeJson(eyeJsonPath, eyeSidecar, true);

files = struct('ieeg', setPath, 'ieegJson', jsonPath, ...
    'channels', channelsPath, 'events', eventsPath, ...
    'electrodes', electrodePath, 'coordinateSystem', coordinatePath, ...
    'eyeTracking', eyePath, 'eyeTrackingJson', eyeJsonPath);
end

function label = normalizeLabel(value, prefix)
label = regexprep(lower(string(value)), '[^a-zA-Z0-9]', '');
label = erase(label, erase(prefix, '-'));
if label == ""
    error('phg:InvalidBidsLabel', 'Empty BIDS label.');
end
label = prefix + label;
end

function value = getOr(input, field, defaultValue)
if isfield(input, field)
    value = input.(field);
else
    value = defaultValue;
end
end

function assertColumns(T, required)
missing = setdiff(required, string(T.Properties.VariableNames));
if ~isempty(missing)
    error('phg:MissingBidsColumns', 'Table lacks: %s', strjoin(missing, ', '));
end
end
