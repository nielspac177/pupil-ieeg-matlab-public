function outputs = runReceptorAlignment(channel, cfg)
%RUNRECEPTORALIGNMENT Does coupling direction track neuromodulator density?
%
%   The paper says it cannot tell whether the polarity dissociation is
%   noradrenergic, cholinergic or neither. This is the cheapest test of that
%   question available without new data: every contact already has an MNI
%   coordinate, and group-average PET maps of the norepinephrine transporter
%   and of the vesicular acetylcholine transporter are public, so each contact
%   can be assigned a normative density and the direction of its coupling
%   regressed on it.
%
%   Three cautions belong with any result from this function.
%
%   1. The PET maps are group averages from other cohorts. No participant here
%      contributed receptor data, so this is a normative-atlas comparison and
%      not a measurement of these brains.
%   2. Receptor maps are spatially smooth and so is the brain, which makes
%      naive p-values badly anticonservative. A variogram-matched spatial null
%      is run for every model and it is the number that counts; the parametric
%      p-values are retained only to show the size of the inflation.
%
%      The dichotomised analysis, on the 285 gated contacts, was negative in
%      the version that matters: NET cleared the spatial null over all
%      contacts (0.044) but not in neocortex alone (0.138), which is where the
%      cortex-versus-subcortex confound has been removed. On that evidence
%      this function used to say a noradrenergic account was not supported.
%
%      On the continuous outcome, which needs no gate and so runs on 913
%      contacts and 750 neocortical ones, the confound-free model does clear
%      the null. The direction is the same and the sample is three times
%      larger. VAChT does not clear it in neocortex alone, so the result is
%      specific to the noradrenergic marker rather than a property of any
%      smooth atlas. Current numbers are written to
%      receptor_coupling_continuous.csv rather than quoted here, so that this
%      comment cannot drift away from them.
%   3. Subcortical PET is degraded by partial-volume effects, and the
%      hippocampus sits next to the ventricles. A hippocampus-versus-cortex
%      contrast on these maps can be an imaging artefact, which is why the
%      neocortex-only model is reported alongside: it removes that confound,
%      and only an effect that survives it is worth discussing.
%
%   Maps are not downloaded automatically. Place them in <workDir>/cache/pet
%   from https://github.com/netneurolab/hansen_receptors (Hansen et al. 2022,
%   Nature Neuroscience); the function reports what is missing and returns
%   empty rather than fetching anything during a pipeline run.

arguments
    channel table
    cfg (1,1) struct
end

petDir = fullfile(cfg.workDir, 'cache', 'pet');
wanted = struct( ...
    'name', {'NET', 'VAChT'}, ...
    'file', {'NAT_MRB_hc77_ding.nii', 'VAChT_feobv_hc18_aghourian_sum.nii'}, ...
    'label', {'Norepinephrine transporter', 'Vesicular ACh transporter'});

present = arrayfun(@(m) isfile(fullfile(petDir, m.file)), wanted);
if ~any(present)
    fprintf(['[PHG] Receptor alignment skipped: no PET maps in %s. ' ...
        'See runReceptorAlignment.m for the source.\n'], petDir);
    outputs = struct('resultTable', table.empty, 'contactTable', table.empty);
    return
end

coordinates = channel.XYZMNI;
frame = table;
frame.PtID = categorical(string(channel.PtID));
frame.Lead = categorical(phg.parseLeadLabel(string(channel.Label)));
frame.Region = phg.cleanRegionLabels(channel.NMM);
frame.HasExcursion = channel.RespAreaNet ~= 0;
frame.IsDilation = double(channel.RespAreaNet > 0);

markerNames = strings(0, 1);
for m = 1:numel(wanted)
    if ~present(m)
        continue
    end
    density = localSampleMap(fullfile(petDir, wanted(m).file), coordinates);
    % z-scored across sampled contacts so that odds ratios are per standard
    % deviation of density and are comparable between markers.
    frame.(wanted(m).name) = (density - mean(density, 'omitnan')) ./ ...
        std(density, 'omitnan');
    markerNames(end+1, 1) = string(wanted(m).name); %#ok<AGROW>
end

excursion = frame(frame.HasExcursion, :);
excursion = excursion(all(isfinite(excursion{:, markerNames}), 2), :);
excursion.PtID = removecats(excursion.PtID);

% Neocortex-only replicate: the mesial temporal structures are dropped, so a
% surviving effect cannot be the hippocampus-versus-cortex contrast in disguise.
neocortex = excursion(~ismember(excursion.Region, ["Hippocampus", "Amygdala"]), :);
neocortex.PtID = removecats(neocortex.PtID);

rows = cell(0, 1);
for k = 1:numel(markerNames)
    marker = markerNames(k);
    rows{end+1, 1} = localFit(excursion, marker, "All excursion contacts"); %#ok<AGROW>
    rows{end+1, 1} = localFit(neocortex, marker, "Neocortex only"); %#ok<AGROW>
end
resultTable = vertcat(rows{:});
phg.writeTableAtomic(resultTable, ...
    fullfile(cfg.tableDir, 'receptor_alignment.csv'));

%% ------------------------------- Continuous outcome, with the spatial null
% The analysis above dichotomises coupling direction and is restricted to the
% contacts that cleared the excursion gate, which is the outcome this paper
% used before the primary moved to a continuous signed amplitude. Repeating it
% on the current primary is the same correction applied to the replication:
% the outcome needs no gate, so the confound-free neocortex-only test runs on
% 750 contacts rather than 222.
%
% The spatial null is computed here rather than in a separate driver script,
% because it is the number that decides the result and a number that decides a
% result should not live outside the pipeline that produces the paper.
scale = median(abs(channel.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
frame.Coupling = asinh(channel.FitHeight ./ scale);
frame.CouplingSd = frame.Coupling ./ std(frame.Coupling, 'omitnan');

continuousRows = cell(0, 1);
for k = 1:numel(markerNames)
    marker = markerNames(k);
    for scope = ["All contacts", "Neocortex only"]
        keep = isfinite(frame.(marker)) & isfinite(frame.Coupling);
        if scope == "Neocortex only"
            keep = keep & ~ismember(frame.Region, ["Hippocampus", "Amygdala"]);
        end
        subset = frame(keep, :);
        subset.PtID = removecats(subset.PtID);
        subset.Value = subset.(marker);
        continuousRows{end + 1, 1} = localContinuousFit(subset, ...
            coordinates(keep, :), marker, scope, ...
            cfg.receptor.numSpatialSurrogates); %#ok<AGROW>
    end
end
continuousTable = vertcat(continuousRows{:});

% Both markers and both scopes form one family of four tests, so the decision
% is made on Benjamini-Hochberg q values rather than on the raw spatial-null
% p-values. Reporting the smallest of four uncorrected p-values as though it
% stood alone is the error this column exists to prevent.
continuousTable.fdr_q_value = phg.benjaminiHochberg( ...
    continuousTable.spatial_null_p_value);
phg.writeTableAtomic(continuousTable, ...
    fullfile(cfg.tableDir, 'receptor_coupling_continuous.csv'));

for k = 1:height(continuousTable)
    fprintf(['[PHG] Receptor, continuous outcome: %s, %s — beta %+.3f ' ...
        '(%.2f SD), parametric P = %.3g, spatial-null P = %.3f, ' ...
        'q = %.3f (n = %d).\n'], continuousTable.marker(k), ...
        continuousTable.sample(k), continuousTable.estimate(k), ...
        continuousTable.standardised(k), ...
        continuousTable.parametric_p_value(k), ...
        continuousTable.spatial_null_p_value(k), ...
        continuousTable.fdr_q_value(k), continuousTable.n_contacts(k));
end

outputs = struct('resultTable', resultTable, 'contactTable', excursion, ...
    'continuousTable', continuousTable);

for k = 1:height(resultTable)
    fprintf(['[PHG] Receptor alignment (parametric P shown; see ' ...
        'receptor_spatial_null.csv for the number that counts): ' ...
        '%s, %s — OR %.2f [%.2f, %.2f], P = %.3g (n = %d).\n'], ...
        resultTable.marker(k), resultTable.sample(k), ...
        resultTable.odds_ratio(k), resultTable.odds_ratio_ci95_low(k), ...
        resultTable.odds_ratio_ci95_high(k), resultTable.p_value(k), ...
        resultTable.n_contacts(k));
end
end

% -------------------------------------------------------------------------
function density = localSampleMap(mapFile, coordinates)
%LOCALSAMPLEMAP Nearest-voxel value of an MNI152 volume at each coordinate.
info = niftiinfo(mapFile);
volume = double(niftiread(mapFile));
% niftiinfo returns the 0-based affine transposed; invert it to go from
% millimetres to voxel indices, then shift to MATLAB's 1-based subscripts.
affine = info.Transform.T';
voxel = round(affine \ [coordinates'; ones(1, size(coordinates, 1))])';
voxel = voxel(:, 1:3) + 1;
inside = all(voxel >= 1, 2) & voxel(:,1) <= info.ImageSize(1) & ...
    voxel(:,2) <= info.ImageSize(2) & voxel(:,3) <= info.ImageSize(3);
density = nan(size(coordinates, 1), 1);
density(inside) = volume(sub2ind(info.ImageSize, voxel(inside,1), ...
    voxel(inside,2), voxel(inside,3)));
end

function row = localFit(frame, marker, sampleLabel)
model = fitglme(frame, ...
    sprintf('IsDilation ~ %s + (1|PtID) + (1|PtID:Lead)', marker), ...
    'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
coefficients = model.Coefficients;
k = find(string(coefficients.Name) == marker, 1);
row = table(marker, string(sampleLabel), height(frame), ...
    exp(coefficients.Estimate(k)), exp(coefficients.Lower(k)), ...
    exp(coefficients.Upper(k)), coefficients.pValue(k), ...
    'VariableNames', {'marker', 'sample', 'n_contacts', 'odds_ratio', ...
    'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'p_value'});
end

% -------------------------------------------------------------------------
function row = localContinuousFit(subset, coordinates, marker, scope, ...
    nSurrogates)
%LOCALCONTINUOUSFIT Coupling on density, with a variogram-matched null.
%   The observed coefficient is compared against models refitted on surrogate
%   density maps that preserve the spatial smoothness of the real map and
%   destroy only its correspondence with anatomy. A parametric p-value on two
%   smooth maps is anticonservative; this is the number that counts.

formula = 'Coupling ~ Value + (1|PtID) + (1|PtID:Lead)';
model = fitlme(subset, formula, 'FitMethod', 'REML');
[~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
index = find(string(stats.Name) == "Value", 1);
observed = stats.Estimate(index);

surrogates = phg.spatialNullSurrogates(subset.Value, coordinates, nSurrogates);
nullBeta = nan(nSurrogates, 1);
surrogateFrame = subset;
for k = 1:nSurrogates
    surrogateFrame.Value = surrogates(:, k);
    try
        nullModel = fitlme(surrogateFrame, formula, 'FitMethod', 'REML');
        nullBeta(k) = nullModel.Coefficients.Estimate( ...
            strcmp(nullModel.Coefficients.Name, 'Value'));
    catch
    end
end
valid = isfinite(nullBeta);
spatialP = (1 + sum(abs(nullBeta(valid)) >= abs(observed))) / (1 + sum(valid));

row = table(marker, scope, height(subset), observed, stats.Lower(index), ...
    stats.Upper(index), observed / std(subset.Coupling, 'omitnan'), ...
    stats.pValue(index), spatialP, sum(valid), 'VariableNames', ...
    {'marker', 'sample', 'n_contacts', 'estimate', 'ci95_low', 'ci95_high', ...
    'standardised', 'parametric_p_value', 'spatial_null_p_value', ...
    'n_surrogates'});
end
