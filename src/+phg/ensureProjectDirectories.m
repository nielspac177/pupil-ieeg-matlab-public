function ensureProjectDirectories(cfg)
%ENSUREPROJECTDIRECTORIES Create reproducible output directories.

dirs = [cfg.resultDir, cfg.figureDir, cfg.tableDir, cfg.workDir, ...
    cfg.cacheDir, cfg.logDir];
for k = 1:numel(dirs)
    if ~isfolder(dirs(k))
        mkdir(dirs(k));
    end
end
end
