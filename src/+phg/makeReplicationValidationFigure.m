function qc = makeReplicationValidationFigure(cfg, options)
%MAKEREPLICATIONVALIDATIONFIGURE The prespecified pre-analysis sanity check.
%
% docs/replication_plan_ebrains.md section 5.2 forbids fitting any model before
% the ingested signals have been looked at: the pupil trace must look like a
% pupil trace, and the high-gamma envelope must have peaks that look like
% peaks. This figure is that gate. It is produced before H1 is estimated and it
% shows raw signal, not a summary, because a summary can look reasonable while
% the underlying trace is nonsense.

arguments
    cfg (1,1) struct
    options.runDir (1,1) string = ""
    options.channel (1,1) string = ""
    options.stem (1,1) string = "FigR1_replication_validation"
end

if options.runDir == ""
    candidates = dir(fullfile(cfg.replication.stageDir, 'sub-*'));
    candidates = candidates([candidates.isdir]);
    options.runDir = string(fullfile(cfg.replication.stageDir, ...
        candidates(1).name));
end

run = phg.loadReplicationRun(options.runDir, 'includeRipple', false);
style = phg.figureStyle;

% Prefer a hippocampal contact: it is the one the primary contrast depends on.
if options.channel == ""
    hippocampal = find(run.channelRegion == "hippocampus" & ...
        run.inAnalysisFrame, 1);
    if isempty(hippocampal)
        hippocampal = find(run.inAnalysisFrame, 1);
    end
else
    hippocampal = find(run.channelName == options.channel, 1);
end

[power, threshold, peakTimes] = localEnvelope(run, hippocampal, cfg);

figureHandle = figure('Units', 'inches', 'Position', [1 1 7.2 6.4], ...
    'Color', 'w');
layout = tiledlayout(figureHandle, 3, 2);
phg.setSafeLayout(layout);

% ---------------------------------------------------------------- panel a
ax = nexttile(layout, [1 2]);
window = localPickWindow(run, peakTimes);
inWindow = run.pupilTime >= window(1) & run.pupilTime <= window(2);
% Both traces are shown about their own means. The comparison this panel has
% to support is one of shape -- does gaze regression distort the pupil signal
% -- and an absolute offset between the two would squash both into ribbons.
rawTrace = run.pupil(inWindow) - mean(run.pupil(inWindow), 'omitnan');
regressedTrace = run.pupilGazeRegressed(inWindow) - ...
    mean(run.pupilGazeRegressed(inWindow), 'omitnan');
hold(ax, 'on');
plot(ax, run.pupilTime(inWindow), rawTrace, ...
    'Color', style.gray, 'LineWidth', style.lineWidth);
plot(ax, run.pupilTime(inWindow), regressedTrace, ...
    'Color', style.constriction, 'LineWidth', style.lineWidth);
events = run.eventOnsets(run.eventOnsets >= window(1) & ...
    run.eventOnsets <= window(2));
yLimits = ylim(ax);
for k = 1:numel(events)
    plot(ax, [events(k) events(k)], yLimits, ':', ...
        'Color', style.lightGray, 'LineWidth', 0.5);
end
ylim(ax, yLimits);
hold(ax, 'off');
xlabel(ax, 'Time in session (s)');
ylabel(ax, 'Pupil area, mean-removed (px^2)');
title(ax, sprintf(['%s: pupil before (grey) and after (dark) gaze ' ...
    'regression; dotted lines are task events'], run.tag), ...
    'Interpreter', 'none');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'a');

% ---------------------------------------------------------------- panel b
ax = nexttile(layout, [1 2]);
inWindow = run.ieegTime >= window(1) & run.ieegTime <= window(2);
hold(ax, 'on');
plot(ax, run.ieegTime(inWindow), power(inWindow), ...
    'Color', style.gray, 'LineWidth', 0.6);
plot(ax, [window(1) window(2)], [threshold threshold], '--', ...
    'Color', style.dilation, 'LineWidth', style.lineWidth);
shown = peakTimes(peakTimes >= window(1) & peakTimes <= window(2));
plot(ax, shown, repmat(threshold, size(shown)), 'v', ...
    'MarkerSize', 3.5, 'MarkerFaceColor', style.dilation, ...
    'MarkerEdgeColor', 'none');
hold(ax, 'off');
xlim(ax, window);
xlabel(ax, 'Time in session (s)');
ylabel(ax, 'High-gamma power (\muV^2)');
% Contact labels contain underscores, which the TeX interpreter would render
% as subscripts (LB_01 becomes LB<sub>0</sub>1). The interpreter is still
% needed for \times and \mu, so the label is escaped rather than disabled.
safeName = strrep(run.channelName(hippocampal), '_', '\_');
title(ax, sprintf(['%s (%s): 70-170 Hz power, threshold at median + 5 ' ...
    '\\times IQR, %d peaks in view'], safeName, ...
    run.channelRegion(hippocampal), numel(shown)), 'Interpreter', 'tex');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'b');

% ---------------------------------------------------------------- panel c
ax = nexttile(layout);
finiteMask = isfinite(run.pupil);
valid = run.pupil(finiteMask);
histogram(ax, valid, 40, 'FaceColor', style.constriction, ...
    'EdgeColor', 'none');
xlabel(ax, 'Pupil area (px^2)');
ylabel(ax, 'Samples');
% Missingness is quoted over the eye-tracker's own coverage. The staged trace
% is padded past the first and last task event, and counting that padding as
% missing data would overstate the loss.
first = find(finiteMask, 1, 'first');
last = find(finiteMask, 1, 'last');
withinCoverage = mean(~finiteMask(first:last));
title(ax, sprintf('Pupil distribution (%.0f%% missing in coverage)', ...
    100 * withinCoverage));
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'c');

% ---------------------------------------------------------------- panel d
ax = nexttile(layout);
counts = localPeakCounts(run, cfg);
histogram(ax, counts, 20, 'FaceColor', style.gray, 'EdgeColor', 'none');
xlabel(ax, 'High-gamma peaks per contact');
ylabel(ax, 'Contacts');
title(ax, sprintf('Peak yield across %d contacts (median %.0f)', ...
    numel(counts), median(counts)));
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'd');

qc = phg.exportPublicationFigure(figureHandle, cfg, options.stem);
close(figureHandle);
end

% -------------------------------------------------------------------------
function [power, threshold, peakTimes] = localEnvelope(run, channel, cfg)
%LOCALENVELOPE Discovery signal chain for one contact.

signalUv = localLaplacian(run, channel);
signalUv(~isfinite(signalUv)) = 0;
nyquist = run.ieegFs / 2;

for f = cfg.replication.notchHz
    if f >= nyquist
        continue
    end
    [b, a] = iirnotch(f / nyquist, cfg.replication.notchBandwidth);
    signalUv = filtfilt(b, a, signalUv);
end
[b, a] = butter(cfg.replication.filterOrder, ...
    cfg.replication.highGammaHz / nyquist, 'bandpass');
power = abs(hilbert(filtfilt(b, a, signalUv))) .^ 2;

inEpoch = power(run.coreMask);
threshold = median(inEpoch) + iqr(inEpoch) * cfg.replication.thresholdIqr;

above = power > threshold & run.coreMask;
edges = diff([false; above; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
peakTimes = nan(numel(starts), 1);
for k = 1:numel(starts)
    [~, offset] = max(power(starts(k):stops(k)));
    peakTimes(k) = run.ieegTime(starts(k) + offset - 1);
end
end

% -------------------------------------------------------------------------
function trace = localLaplacian(run, channel)
sameShaft = find(run.channelShaft == run.channelShaft(channel));
neighbours = sameShaft( ...
    abs(run.channelIndex(sameShaft) - run.channelIndex(channel)) == 1);
if isempty(neighbours)
    trace = run.ieeg(:, channel);
else
    trace = run.ieeg(:, channel) - mean(run.ieeg(:, neighbours), 2);
end
end

% -------------------------------------------------------------------------
function window = localPickWindow(run, peakTimes)
%LOCALPICKWINDOW A 60 s stretch of eligible recording containing peaks.

span = 60;
if isempty(peakTimes)
    window = [run.ieegTime(1), run.ieegTime(1) + span];
    return
end
centre = median(peakTimes);
window = [centre - span / 2, centre + span / 2];
window(1) = max(window(1), run.ieegTime(1));
window(2) = min(window(1) + span, run.ieegTime(end));
end

% -------------------------------------------------------------------------
function counts = localPeakCounts(run, cfg)
%LOCALPEAKCOUNTS Peak yield for every analysis contact in this run.

index = find(run.inAnalysisFrame);
counts = nan(numel(index), 1);
for k = 1:numel(index)
    [~, ~, times] = localEnvelope(run, index(k), cfg);
    counts(k) = numel(times);
end
end
