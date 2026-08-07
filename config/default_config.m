function cfg = default_config(repoRoot)
%DEFAULT_CONFIG Versioned defaults for the pupil-iEEG pipeline.

arguments
    repoRoot (1,1) string = string(fileparts(fileparts(mfilename('fullpath'))))
end

cfg = struct;
cfg.repoRoot = repoRoot;
cfg.legacyRoot = string(fileparts(repoRoot));
cfg.derivedMatFile = fullfile(cfg.legacyRoot, 'PupilHG', 'PHG_20s.mat');
cfg.legacyCodeDir = fullfile(cfg.legacyRoot, 'PupilHG');

cfg.resultDir = fullfile(repoRoot, 'results');
cfg.figureDir = fullfile(cfg.resultDir, 'figures');
cfg.tableDir = fullfile(cfg.resultDir, 'tables');
cfg.workDir = fullfile(repoRoot, 'work');
cfg.cacheDir = fullfile(cfg.workDir, 'cache');
cfg.logDir = fullfile(cfg.workDir, 'logs');
cfg.bidsRoot = fullfile(cfg.workDir, 'bids');

cfg.randomSeed = 20260802;
cfg.legacySelectionThreshold = 0.10;
cfg.topRegionCount = 10;

% Transparent-brain rendering. Point this to LeGUI/SPM TPM.nii locally.
cfg.leGUIRoot = "";
cfg.templateTpmFile = "";
cfg.templateSurfaceCache = fullfile(cfg.cacheDir, 'mni_transparent_surface.mat');
cfg.brainAlpha = 0.16;
cfg.allContactColor = [0.58 0.60 0.63];
cfg.allContactAlpha = 0.30;
cfg.selectedContactSize = 26;
cfg.allContactSize = 7;

% BIDS staging. The repository never publishes this directory automatically.
cfg.bidsVersion = "1.11.1";
cfg.createBidsScaffold = true;
cfg.rawManifestFile = fullfile(repoRoot, 'config', 'raw_manifest.tsv');
cfg.datasetName = "Pupil-linked human intracranial electrophysiology";

% Ripple detector defaults. The event detector operates on all eligible
% hippocampal contacts; it must never be restricted by legacy pupil selection.
cfg.ripple = struct;
cfg.ripple.bandHz = [80 120];
cfg.ripple.filterOrder = 4;
cfg.ripple.rmsWindowSeconds = 0.020;
cfg.ripple.minimumZ = 2.5;
cfg.ripple.maximumZ = 9.0;
cfg.ripple.minimumDurationSeconds = 0.038;
cfg.ripple.maximumRawAmplitudeUv = 300;
cfg.ripple.spectralWindowSeconds = 0.250;

% Pupil prediction defaults based on the EEGNet regression precedent.
cfg.prediction = struct;
cfg.prediction.trainingFraction = 0.80;
cfg.prediction.initialLearnRate = 1e-5;
cfg.prediction.maxEpochs = 10000;
cfg.prediction.validationPatience = 500;
cfg.prediction.minimumDelta = 1e-3;
cfg.prediction.temporalFilters = 8;
cfg.prediction.depthMultiplier = 2;
cfg.prediction.dropout = 0.50;
end
