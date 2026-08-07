%% RUN_DEMO Everything in this repository that runs without participant data.
% The full pipeline (run_all) needs the derived dataset, which is not public.
% This script exercises the parts that do not: the methods schematic, which is
% drawn from synthetic traces, and the ripple detector, on a simulated signal.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'config'));

cfg = default_config(repoRoot);
phg.ensureProjectDirectories(cfg);

fprintf('[PHG] Methods schematic (synthetic traces)\n');
phg.makeMethodsFigure(cfg);

fprintf('[PHG] Ripple detector demonstration (simulated signal)\n');
run(fullfile(repoRoot, 'examples', 'demo_ripple_detector.m'));

fprintf('[PHG] Demo complete. Output in %s\n', cfg.figureDir);
