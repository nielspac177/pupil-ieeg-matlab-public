%% RUN_ALL Reproduce all currently available analyses and figures.
% Raw-dependent stages are reported as skipped until their required inputs are
% configured. No Python or R code is called by this pipeline.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'config'));

cfg = default_config(repoRoot);
localConfig = fullfile(repoRoot, 'config', 'local_config.m');
if isfile(localConfig)
    run(localConfig);
end

rng(cfg.randomSeed, 'twister');
phg.ensureProjectDirectories(cfg);

fprintf('[PHG] MATLAB-only pipeline started: %s\n', char(string(datetime('now'))));
fprintf('[PHG] Repository: %s\n', char(repoRoot));

[tables, audit] = phg.loadDerivedTables(cfg);
phg.writeTableAtomic(audit, fullfile(cfg.tableDir, 'derived_data_audit.csv'));

phg.makeMethodsFigure(cfg);

analysis = phg.runDerivedAnalyses(tables.channel, cfg);
ripple = phg.runRippleBandTest(tables.wavelet, tables.channel, cfg);
receptors = phg.runReceptorAlignment(tables.channel, cfg);

phg.makeDerivedPublicationFigures(tables, cfg, analysis);
phg.makeElectrodeLocalization(tables.channel, cfg, 'mode', "effect");
phg.makeElectrodeLocalization(tables.channel, cfg, 'mode', "region");
phg.makeElectrodeLocalization(tables.channel, cfg, 'mode', "polarity");

rawStatus = phg.auditRawAvailability(cfg);
phg.writeTableAtomic(rawStatus, fullfile(cfg.tableDir, 'raw_data_readiness.csv'));

if cfg.createBidsScaffold
    phg.createBidsScaffold(cfg);
end

phg.writeRunManifest(cfg, rawStatus);
fprintf('[PHG] Pipeline completed: %s\n', char(string(datetime('now'))));
