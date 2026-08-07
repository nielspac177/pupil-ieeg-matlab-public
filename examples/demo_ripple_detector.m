%% DEMO_RIPPLE_DETECTOR Verify the detector on a known synthetic burst.
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

fs = 1000;
time = (0:fs*4-1)' / fs;
rng(7, 'twister');
signal = 8 * randn(size(time));
burst = time >= 1.95 & time < 2.05;
signal(burst) = signal(burst) + 75 * sin(2*pi*100*time(burst));

events = phg.detectRipples(signal, fs, ...
    'minimumZ', 2.0, 'maximumZ', 20, 'maximumRawAmplitudeUv', 300);
accepted = events(events.accepted, :);
assert(~isempty(accepted), 'Synthetic 100-Hz burst was not detected.');
assert(any(abs(accepted.peak_frequency_hz - 100) < 8), ...
    'Detected event did not peak near 100 Hz.');
disp(accepted);
