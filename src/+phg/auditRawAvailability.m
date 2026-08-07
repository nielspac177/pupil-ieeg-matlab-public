function status = auditRawAvailability(cfg)
%AUDITRAWAVAILABILITY Report raw-data readiness without exporting private paths.

if ~isfile(cfg.rawManifestFile)
    error('phg:MissingRawManifest', ...
        'Raw-data manifest not found: %s', cfg.rawManifestFile);
end

manifest = readtable(cfg.rawManifestFile, 'FileType', 'text', ...
    'Delimiter', '\t', 'TextType', 'string');
pathColumns = ["raw_ieeg", "nev", "pupil", "channel_map", ...
    "bad_channels", "anatomy", "soz_ied"];
missingColumns = setdiff(pathColumns, string(manifest.Properties.VariableNames));
if ~isempty(missingColumns)
    error('phg:InvalidRawManifest', 'Missing columns: %s', ...
        strjoin(missingColumns, ', '));
end

present = false(height(manifest), numel(pathColumns));
for c = 1:numel(pathColumns)
    value = manifest.(pathColumns(c));
    for r = 1:height(manifest)
        if ~ismissing(value(r)) && strlength(value(r)) > 0
            present(r,c) = isfile(value(r)) || isfolder(value(r));
        end
    end
end

coreIndex = ismember(pathColumns, ["raw_ieeg", "pupil", "channel_map"]);
confirmatoryIndex = ismember(pathColumns, ...
    ["raw_ieeg", "pupil", "channel_map", "bad_channels", "soz_ied"]);
bidsIndex = ismember(pathColumns, ...
    ["raw_ieeg", "pupil", "channel_map", "anatomy"]);

status = manifest(:, {'participant_id', 'session_id', 'task'});
status.n_files_present = sum(present, 2);
status.core_analysis_ready = all(present(:, coreIndex), 2);
status.confirmatory_analysis_ready = all(present(:, confirmatoryIndex), 2);
status.bids_release_candidate = all(present(:, bidsIndex), 2);
status.supplied_archive_status = repmat( ...
    "raw_inputs_not_available", height(status), 1);
status.supplied_archive_status(status.core_analysis_ready) = ...
    "core_inputs_available";

fprintf('[PHG] Raw readiness: %d/%d sessions have core inputs.\n', ...
    sum(status.core_analysis_ready), height(status));
if ~any(status.core_analysis_ready)
    fprintf(['[PHG] Conclusion: raw-dependent analyses are not estimable ' ...
        'from the supplied archive.\n']);
end
end
