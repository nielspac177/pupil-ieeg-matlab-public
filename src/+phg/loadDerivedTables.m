function [tables, audit] = loadDerivedTables(cfg)
%LOADDERIVEDTABLES Load and validate the authoritative legacy MAT object.

arguments
    cfg (1,1) struct
end

if ~isfile(cfg.derivedMatFile)
    error('phg:MissingDerivedData', ...
        'Derived MAT file not found: %s', cfg.derivedMatFile);
end
if isfolder(cfg.legacyCodeDir)
    addpath(cfg.legacyCodeDir);
end

loaded = load(cfg.derivedMatFile);
if ~isfield(loaded, 'obj')
    error('phg:InvalidDerivedData', ...
        'Expected an object named obj in %s.', cfg.derivedMatFile);
end

legacy = loaded.obj;
tables = struct;
tables.channel = legacy.MultPtTable;
tables.channelPair = legacy.MultPtTableHG;
tables.wavelet = legacy.MultPtTableWav;

required = ["PtID", "Chan", "NMM", "XYZMNI", "RespSig", ...
    "RespAreaNet", "RespAreaAbs", "FitShift", "FitR2"];
missing = setdiff(required, string(tables.channel.Properties.VariableNames));
if ~isempty(missing)
    error('phg:InvalidDerivedData', ...
        'Channel table lacks required variables: %s', strjoin(missing, ', '));
end

T = tables.channel;
isLegacySelected = T.RespSig > cfg.legacySelectionThreshold;
% Match the hippocampus proper. A substring match on 'Hipp' also captures
% 'PHG parahippocampal gyrus', which is a distinct parcel.
isHippocampus = phg.cleanRegionLabels(T.NMM) == "Hippocampus";
hasCoordinates = ~any(isnan(T.XYZMNI), 2);

metric = [
    "patients"
    "channels"
    "legacy_selected_channels"
    "legacy_selection_threshold"
    "channels_with_mni_coordinates"
    "hippocampal_channels"
    "patients_with_hippocampal_coverage"
    "legacy_selected_hippocampal_channels"
    "channel_pair_rows"
    "wavelet_rows"
    ];
value = [
    numel(unique(string(T.PtID)))
    height(T)
    sum(isLegacySelected)
    cfg.legacySelectionThreshold
    sum(hasCoordinates)
    sum(isHippocampus)
    numel(unique(string(T.PtID(isHippocampus))))
    sum(isHippocampus & isLegacySelected)
    height(tables.channelPair)
    height(tables.wavelet)
    ];
audit = table(metric, value, 'VariableNames', {'metric', 'value'});

fprintf('[PHG] Loaded %d channels from %d patients; %d legacy-selected.\n', ...
    height(T), numel(unique(string(T.PtID))), sum(isLegacySelected));
end
