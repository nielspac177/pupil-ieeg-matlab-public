function run = loadReplicationRun(runDir, options)
%LOADREPLICATIONRUN Read one staged EBRAINS run into MATLAB.
%   The Python ingestion step (tools/ingest_ebrains.py) decodes MEF3, fixes the
%   shared time base and decimates to 1000 Hz. It performs no measurement. This
%   function reads what it wrote and hands back arrays in MATLAB orientation.
%
%   HDF5 is row-major and MATLAB is column-major, so a samples-by-channels
%   array written by h5py arrives as channels-by-samples. Every matrix is
%   transposed here, once, so that no caller has to remember it.

arguments
    runDir (1,1) string
    % The 4 kHz hippocampal arrays are roughly a third of a run's bulk and are
    % only wanted by the ripple analysis. Everything else should decline them:
    % holding them made the full-cohort measurement exhaust memory.
    options.includeRipple (1,1) logical = true
end

h5File = fullfile(runDir, 'run.h5');
metaFile = fullfile(runDir, 'meta.json');
if ~isfile(h5File) || ~isfile(metaFile)
    error('phg:MissingReplicationRun', ...
        'Staged run is incomplete: %s', runDir);
end

meta = jsondecode(fileread(metaFile));

run = struct;
run.tag = string(meta.tag);
run.subject = string(meta.subject);
run.task = string(meta.task);
run.runIndex = string(meta.run);
run.ieegFs = meta.ieeg_fs;
run.pupilFs = meta.pupil_fs;
run.anchorValid = logical(meta.anchor_valid);

% Channel metadata. jsondecode collapses a struct array into a struct of
% arrays only when every entry has identical fields, which it does here.
channels = meta.channels;
run.channelName = string({channels.name})';
run.channelRegion = string({channels.region})';
run.channelShaft = string({channels.shaft})';
run.channelIndex = cellfun(@double, {channels.contact_index})';
run.inAnalysisFrame = logical(cellfun(@double, {channels.in_analysis_frame})');
run.channelXyz = [ ...
    cellfun(@(v) str2double(string(v)), {channels.x})', ...
    cellfun(@(v) str2double(string(v)), {channels.y})', ...
    cellfun(@(v) str2double(string(v)), {channels.z})'];

run.ieeg = h5read(h5File, '/ieeg')';                 % samples x channels
run.ieegTime = double(h5read(h5File, '/ieeg_time'));
run.coreMask = h5read(h5File, '/core_mask') > 0;
run.maskRuleB = h5read(h5File, '/mask_rule_b') > 0;
run.maskRuleC = h5read(h5File, '/mask_rule_c') > 0;
run.pupil = double(h5read(h5File, '/pupil'));
run.pupilGazeRegressed = double(h5read(h5File, '/pupil_gaze_regressed'));
run.pupilTime = double(h5read(h5File, '/pupil_time'));
run.eventOnsets = double(h5read(h5File, '/event_onsets'));

info = h5info(h5File);
staged = string({info.Datasets.Name});
% The channel list is always reported, even when the signal is not loaded, so a
% caller can tell that a run has hippocampal coverage without paying for it.
run.rippleChannels = strings(0, 1);
if isfield(meta, 'ripple_channels') && ~isempty(meta.ripple_channels)
    run.rippleChannels = string(meta.ripple_channels);
    run.rippleChannels(run.rippleChannels == "") = [];
end

run.rippleIeeg = [];
run.rippleTime = [];
run.rippleFs = NaN;
if ismember("ripple_ieeg", staged) && options.includeRipple
    run.rippleIeeg = h5read(h5File, '/ripple_ieeg')';
    run.rippleTime = double(h5read(h5File, '/ripple_time'));
    run.rippleFs = meta.ripple_fs;
end
end
