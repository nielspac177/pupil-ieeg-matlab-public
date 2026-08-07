function writeRunManifest(cfg, rawStatus)
%WRITERUNMANIFEST Record deterministic software and input metadata.

toolboxes = ver;
toolboxInfo = repmat(struct('name', '', 'version', ''), numel(toolboxes), 1);
for k = 1:numel(toolboxes)
    toolboxInfo(k).name = toolboxes(k).Name;
    toolboxInfo(k).version = toolboxes(k).Version;
end

manifest = struct;
manifest.generated_at = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
manifest.matlab_version = version;
manifest.random_seed = cfg.randomSeed;
% The file name and its checksum identify the input; the directory it happened
% to sit in on one machine does not, and this manifest is committed.
[~, derivedName, derivedExt] = fileparts(char(cfg.derivedMatFile));
manifest.derived_mat_file = [derivedName derivedExt];
manifest.derived_mat_sha256 = phg.sha256File(cfg.derivedMatFile);
manifest.raw_sessions_core_ready = sum(rawStatus.core_analysis_ready);
manifest.raw_sessions_total = height(rawStatus);
manifest.toolboxes = toolboxInfo;
phg.writeJson(fullfile(cfg.resultDir, 'run_manifest.json'), manifest, true);
end
