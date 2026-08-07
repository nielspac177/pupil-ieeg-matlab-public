%% RUN_TESTS Execute the MATLAB unit and integration tests.
repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'config'));
suite = testsuite(fullfile(repoRoot, 'tests'), 'IncludeSubfolders', true);
results = run(suite);
disp(table(results));
assert(all([results.Passed]), 'phg:TestsFailed', ...
    'At least one pupil-iEEG test failed.');
