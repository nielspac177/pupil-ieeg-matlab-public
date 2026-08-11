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

% Spatial-gradient analysis. A shaft must span this much distance from the
% hippocampal centroid before within-shaft variation means anything.
cfg.gradient = struct;
cfg.gradient.minimumShaftSpanMm = 5;
cfg.gradient.numPermutations = 300;
% Candidate origins for the gradient-origin test need enough contacts to
% define a centroid, and models within this AIC of the best are treated as
% indistinguishable, following the conventional threshold of about four.
cfg.gradient.minimumOriginContacts = 10;
cfg.gradient.aicIndistinguishable = 4;

% Architecture analyses: frequency profile, response latency, spatial scale.
cfg.architecture = struct;
% Surrogate count for the variogram-matched spatial null used by the receptor
% analysis. This is the number that decides that result, so it lives here
% rather than in a driver script.
cfg.receptor.numSpatialSurrogates = 500;

cfg.architecture.minimumFitR2 = 0.25;
cfg.architecture.maximumLagMs = 20000;
cfg.architecture.distanceBinsMm = [0 5 10 15 20 30 40 60 80 120];

% EBRAINS replication cohort. Every constant here is copied from the discovery
% analysis object rather than chosen, so that a difference in result between
% cohorts is attributable to the data. Verified against PHG_20s.mat:
% obj.Fs = 1000, obj.BandPass = [70 170], obj.RRType = 'Laplace', obj.ThrX = 5,
% obj.TimeRng = [-20000 20000]. See docs/replication_amendment_01.md.
cfg.replication = struct;
% The staged replication cohort is tens of gigabytes and lives on external
% storage. Point cfg.replication.stageDir at it in local_config.m; the default
% below keeps the repository free of machine-local paths.
cfg.replication.stageDir = fullfile(cfg.workDir, 'ebrains_staging');
cfg.replication.ieegFs = 1000;
cfg.replication.pupilFs = 150;
cfg.replication.highGammaHz = [70 170];
cfg.replication.filterOrder = 1;          % butter(1, ...) then filtfilt
cfg.replication.notchHz = 60:60:240;
cfg.replication.notchBandwidth = 0.012;
cfg.replication.thresholdIqr = 5;         % median + 5 * IQR
cfg.replication.timeRangeMs = [-20000 20000];
cfg.replication.responseWindowSeconds = 5;
cfg.replication.numPermutations = 100;
cfg.replication.maxTrialsPerPermutation = 1000;
% Prespecified in the amendment, not inherited from discovery.
cfg.replication.maxWindowMissingFraction = 0.40;
cfg.replication.minPeaksPerContact = 10;
cfg.replication.numCircularShifts = 1000;
% The decision rule, fixed in docs/replication_plan_ebrains.md before the data
% arrived. discoveryInterval is the discovery 95% interval for the primary
% odds ratio and defines the "directionally consistent but underpowered" band.
cfg.replication.discoveryInterval = [0.021 0.207];
cfg.replication.minHippocampalExcursion = 10;
cfg.replication.overlapAdjustmentGap = 0.05;

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
