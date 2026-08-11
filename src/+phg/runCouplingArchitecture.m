function outputs = runCouplingArchitecture(tables, cfg)
%RUNCOUPLINGARCHITECTURE Frequency, timing and spatial structure of coupling.
%
% Three questions the binary primary could not ask, all on the continuous
% signed coupling amplitude.
%
% 1. Frequency. Pfeffer et al. (2022) report in source-resolved MEG that the
%    sign of pupil-activity coupling depends on frequency band, being negative
%    at low frequencies across most of the brain with peaks in anterior
%    hippocampus. If that is right, the hippocampal-versus-neocortical
%    dissociation reported here should be frequency-dependent rather than
%    fixed. The wavelet table carries the same peri-peak measurement at 71
%    frequencies between 2 and 256 Hz, so the prediction is directly testable.
%
% 2. Timing. FitShift is the lag between the high-gamma peak and the fitted
%    pupil response. A single common driver predicts a shared lag; separate
%    pathways do not. This is a null in these data and is reported as one.
%
% 3. Spatial scale. The gradient analysis shows coupling varies with position.
%    That is only interpretable if coupling is spatially organised at all, so
%    the autocorrelation of coupling against inter-contact distance is
%    estimated and a length constant fitted. It also justifies the
%    within-shaft permutation null used for the gradient, by measuring the
%    autocorrelation that null is designed to preserve.
%
% The frequency analysis inherits a limitation the others do not: the wavelet
% table exists only for contacts that passed the historical selection rule, so
% it covers 177 contacts rather than 913. That is stated wherever it is used.

arguments
    tables (1,1) struct
    cfg (1,1) struct
end

T = tables.channel;
W = tables.wavelet;
warningState = warning('off', 'all');

%% ------------------------------------------------ 1. Frequency profile
regionW = phg.cleanRegionLabels(W.NMM);
classW = repmat("Extrahippocampal", height(W), 1);
classW(regionW == "Hippocampus") = "Hippocampal";
scaleW = median(abs(W.FitHeight), 'omitnan');

frequencyFrame = table(categorical(string(W.PtID)), ...
    categorical(string(W.PtID) + "/" + phg.parseLeadLabel(W.Label)), ...
    categorical(string(W.PtID) + "/" + string(W.Label)), ...
    asinh(W.FitHeight ./ scaleW), log2(W.freq), W.freq, ...
    categorical(classW, ["Extrahippocampal", "Hippocampal"]), ...
    'VariableNames', {'PtID', 'Lead', 'Contact', 'Coupling', 'LogFreq', ...
    'Frequency', 'ContactClass'});
frequencyFrame = frequencyFrame(isfinite(frequencyFrame.Coupling), :);

% Every contact contributes one row per wavelet band, so a contact appears 71
% times and those 71 rows come from a single peri-peak curve. Grouping only by
% patient and shaft leaves that repetition unmodelled, which treats correlated
% measurements as independent replicates: the earlier specification returned
% roughly 12,500 denominator degrees of freedom from 177 contacts, and a
% correspondingly extreme p-value.
%
% The term of interest is an interaction between contact class, which varies
% between contacts, and frequency, which varies within them. Testing it
% requires knowing how much the frequency slope varies from contact to
% contact, so the slope is allowed to vary by contact. Satterthwaite then
% returns degrees of freedom on the order of the number of contacts, which is
% the number of independent slopes the contrast is actually estimated from.
interactionModel = fitlme(frequencyFrame, ...
    ['Coupling ~ ContactClass*LogFreq + (1|PtID) + (1|PtID:Lead) + ' ...
     '(1 + LogFreq|Contact)'], 'FitMethod', 'REML');
[~, ~, interactionStats] = fixedEffects(interactionModel, ...
    'DFMethod', 'satterthwaite');

% The superseded specification is fitted and exported alongside, so the size
% of the correction is on the record rather than only its result.
naiveModel = fitlme(frequencyFrame, ...
    'Coupling ~ ContactClass*LogFreq + (1|PtID) + (1|PtID:Lead)', ...
    'FitMethod', 'REML');
[~, ~, naiveStats] = fixedEffects(naiveModel, 'DFMethod', 'satterthwaite');
naiveIndex = find(contains(string(naiveStats.Name), ':'), 1);
interactionIndex = find(contains(string(interactionStats.Name), ':'), 1);
frequencySpecification = table( ...
    ["contact_random_slope"; "no_contact_grouping"], ...
    [interactionStats.Estimate(interactionIndex); ...
     naiveStats.Estimate(naiveIndex)], ...
    [interactionStats.SE(interactionIndex); naiveStats.SE(naiveIndex)], ...
    [interactionStats.DF(interactionIndex); naiveStats.DF(naiveIndex)], ...
    [interactionStats.pValue(interactionIndex); ...
     naiveStats.pValue(naiveIndex)], ...
    'VariableNames', {'specification', 'estimate', 'standard_error', ...
    'satterthwaite_df', 'p_value'});
phg.writeTableAtomic(frequencySpecification, ...
    fullfile(cfg.tableDir, 'coupling_frequency_specification.csv'));
fprintf(['[PHG]   frequency interaction: beta %+.4f, P = %.3g with slopes ' ...
    'varying by contact (df %.0f); P = %.3g without (df %.0f)\n'], ...
    interactionStats.Estimate(interactionIndex), ...
    interactionStats.pValue(interactionIndex), ...
    interactionStats.DF(interactionIndex), ...
    naiveStats.pValue(naiveIndex), naiveStats.DF(naiveIndex));
frequencyModel = table(string(interactionStats.Name), ...
    interactionStats.Estimate, interactionStats.SE, ...
    interactionStats.tStat, interactionStats.DF, interactionStats.pValue, ...
    interactionStats.Lower, interactionStats.Upper, 'VariableNames', ...
    {'term', 'estimate', 'standard_error', 't_statistic', ...
    'satterthwaite_df', 'p_value', 'ci95_low', 'ci95_high'});
phg.writeTableAtomic(frequencyModel, ...
    fullfile(cfg.tableDir, 'coupling_frequency_model.csv'));

edges = [2 4 8 16 32 64 128 256];
bandLabel = strings(0, 1); bandLow = []; bandHigh = [];
hippMedian = []; otherMedian = []; hippN = []; otherN = [];
for b = 1:numel(edges) - 1
    inBand = frequencyFrame.Frequency >= edges(b) & ...
        frequencyFrame.Frequency < edges(b + 1);
    h = inBand & frequencyFrame.ContactClass == "Hippocampal";
    o = inBand & frequencyFrame.ContactClass == "Extrahippocampal";
    if sum(h) < 5 || sum(o) < 5
        continue
    end
    bandLabel(end + 1, 1) = sprintf('%g-%g Hz', edges(b), edges(b + 1));
    bandLow(end + 1, 1) = edges(b);
    bandHigh(end + 1, 1) = edges(b + 1);
    hippMedian(end + 1, 1) = median(frequencyFrame.Coupling(h));
    otherMedian(end + 1, 1) = median(frequencyFrame.Coupling(o));
    hippN(end + 1, 1) = sum(h);
    otherN(end + 1, 1) = sum(o);
end
frequencyBands = table(bandLabel, bandLow, bandHigh, hippN, hippMedian, ...
    otherN, otherMedian, 'VariableNames', {'band', 'low_hz', 'high_hz', ...
    'n_hippocampal', 'hippocampal_median', 'n_extrahippocampal', ...
    'extrahippocampal_median'});
phg.writeTableAtomic(frequencyBands, ...
    fullfile(cfg.tableDir, 'coupling_frequency_bands.csv'));

interaction = frequencyModel(contains(frequencyModel.term, ":"), :);
fprintf(['[PHG] Frequency: class-by-frequency interaction beta = %+.4f, ', ...
    'P = %.3g, on %d contacts\n'], interaction.estimate(1), ...
    interaction.p_value(1), ...
    numel(unique(string(W.PtID) + "|" + string(W.Label))));

%% ------------------------------------------------------- 2. Latency
regionT = phg.cleanRegionLabels(T.NMM);
classT = repmat("Extrahippocampal", height(T), 1);
classT(regionT == "Hippocampus") = "Hippocampal";
usableLag = isfinite(T.FitShift) & T.FitR2 > cfg.architecture.minimumFitR2 & ...
    abs(T.FitShift) < cfg.architecture.maximumLagMs;
lagFrame = table(categorical(string(T.PtID)), ...
    categorical(string(T.PtID) + "/" + phg.parseLeadLabel(T.Label)), ...
    T.FitShift / 1000, ...
    categorical(classT, ["Extrahippocampal", "Hippocampal"]), ...
    'VariableNames', {'PtID', 'Lead', 'LagSeconds', 'ContactClass'});
lagFrame = lagFrame(usableLag, :);

lagModel = fitlme(lagFrame, ...
    'LagSeconds ~ ContactClass + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
[~, ~, lagStats] = fixedEffects(lagModel, 'DFMethod', 'satterthwaite');
index = find(startsWith(string(lagStats.Name), "ContactClass"), 1);
latency = table(height(lagFrame), ...
    median(lagFrame.LagSeconds(lagFrame.ContactClass == "Hippocampal")), ...
    median(lagFrame.LagSeconds(lagFrame.ContactClass == "Extrahippocampal")), ...
    lagStats.Estimate(index), lagStats.Lower(index), lagStats.Upper(index), ...
    lagStats.pValue(index), 'VariableNames', {'n_contacts', ...
    'hippocampal_median_seconds', 'extrahippocampal_median_seconds', ...
    'difference_seconds', 'ci95_low', 'ci95_high', 'p_value'});
phg.writeTableAtomic(latency, ...
    fullfile(cfg.tableDir, 'coupling_latency.csv'));
fprintf(['[PHG] Latency: hippocampal %+.2f s vs extrahippocampal %+.2f s, ', ...
    'difference %+.2f s, P = %.3g\n'], latency.hippocampal_median_seconds, ...
    latency.extrahippocampal_median_seconds, latency.difference_seconds, ...
    latency.p_value);

%% ------------------------------------------------- 3. Spatial scale
coordinates = T.XYZMNI;
scaleT = median(abs(T.FitHeight), 'omitnan');
coupling = asinh(T.FitHeight ./ scaleT);
good = all(isfinite(coordinates), 2) & isfinite(coupling);
patient = string(T.PtID);

binEdges = cfg.architecture.distanceBinsMm;
binCentre = []; binPairs = []; binCorrelation = [];
for b = 1:numel(binEdges) - 1
    left = []; right = [];
    for p = unique(patient(good))'
        index = find(good & patient == p);
        if numel(index) < 3
            continue
        end
        separation = squareform(pdist(coordinates(index, :)));
        [i, j] = find(triu(separation >= binEdges(b) & ...
            separation < binEdges(b + 1), 1));
        left = [left; coupling(index(i))]; %#ok<AGROW>
        right = [right; coupling(index(j))]; %#ok<AGROW>
    end
    if numel(left) < 20
        continue
    end
    binCentre(end + 1, 1) = mean(binEdges(b:b + 1)); %#ok<AGROW>
    binPairs(end + 1, 1) = numel(left); %#ok<AGROW>
    binCorrelation(end + 1, 1) = corr(left, right); %#ok<AGROW>
end

% Exponential decay gives a length constant: the distance over which coupling
% similarity falls by 1/e.
lengthConstant = NaN;
positive = binCorrelation > 0;
if sum(positive) >= 3
    fit = polyfit(binCentre(positive), log(binCorrelation(positive)), 1);
    lengthConstant = -1 / fit(1);
end
spatialScale = table(binCentre, binPairs, binCorrelation, ...
    repmat(lengthConstant, numel(binCentre), 1), 'VariableNames', ...
    {'distance_mm', 'n_pairs', 'correlation', 'length_constant_mm'});
phg.writeTableAtomic(spatialScale, ...
    fullfile(cfg.tableDir, 'coupling_spatial_scale.csv'));
fprintf(['[PHG] Spatial scale: correlation %.2f below %g mm, %.2f at ', ...
    '%g mm; length constant %.1f mm\n'], binCorrelation(1), binEdges(2), ...
    binCorrelation(end), binCentre(end), lengthConstant);

warning(warningState);
outputs = struct('frequencyModel', frequencyModel, ...
    'frequencyBands', frequencyBands, 'latency', latency, ...
    'spatialScale', spatialScale);
end
