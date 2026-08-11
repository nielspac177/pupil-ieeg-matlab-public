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

% Primary analysis, on the scale this field uses: a continuous signed coupling
% value modelled with linear mixed effects, no binarisation and no selection.
% See runContinuousCouplingAnalyses for why this replaced the binary primary.
continuous = phg.runContinuousCouplingAnalyses(tables.channel, cfg);
gradientResult = phg.runSpatialGradient(tables.channel, cfg);
phg.makeContinuousPrimaryFigure(tables.channel, cfg, ...
    'wavelet', tables.wavelet);
originSpecificity = phg.runGradientOriginSpecificity(tables.channel, cfg);
continuousRobustness = phg.runContinuousRobustness(tables.channel, cfg);
effectMeasures = phg.runEffectMeasureComparison(tables.channel, cfg);
architecture = phg.runCouplingArchitecture(tables, cfg);

% The replication judged on the current primary scale, and the cross-cohort
% comparison that decides whether its null is informative or merely quiet.
% Both depend only on stored per-contact measures, so they run here rather
% than inside the replication stage, which needs the raw recordings.
if isfile(fullfile(cfg.tableDir, ...
        'replication_contact_measures_gaze_regressed_ruleA.csv'))
    replicationContinuous = phg.runReplicationContinuous(cfg);
    comparability = phg.runReplicationComparability(tables.channel, cfg);
else
    fprintf(['[PHG] Replication measures absent; skipping the continuous ' ...
        'replication and the cross-cohort comparison.\n']);
end

% The hurdle decomposition is retained as a secondary, comparability analysis:
% it is what earlier versions reported, and the two framings should be visible
% side by side rather than one silently replacing the other.
analysis = phg.runDerivedAnalyses(tables.channel, cfg);
sensitivity = phg.runPrimarySensitivity(tables.channel, cfg);
referenceChecks = phg.runReferenceClassChecks(tables.channel, cfg);
ripple = phg.runRippleBandTest(tables.wavelet, tables.channel, cfg);
receptors = phg.runReceptorAlignment(tables.channel, cfg);
phg.makeReceptorSurfaceFigure(tables.channel, cfg);

phg.makeDerivedPublicationFigures(tables, cfg, analysis);
phg.makeForestTable(cfg, analysis);
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
