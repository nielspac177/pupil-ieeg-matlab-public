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
%      naive p-values anticonservative. The patient and shaft random effects
%      absorb some of that, but a spatially constrained null (spin test or
%      variogram-matched surrogates) is the proper control and is NOT done
%      here. Treat p-values as descriptive.
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

outputs = struct('resultTable', resultTable, 'contactTable', excursion);

for k = 1:height(resultTable)
    fprintf(['[PHG] Receptor alignment (exploratory, no spatial null): ' ...
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
