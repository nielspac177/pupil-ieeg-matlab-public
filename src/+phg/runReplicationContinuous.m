function outputs = runReplicationContinuous(cfg, options)
%RUNREPLICATIONCONTINUOUS Judge the replication on the current primary scale.
%
% The replication was run and reported when the primary outcome was a binary
% direction, and was inconclusive on that outcome. The primary has since become
% a continuous signed coupling amplitude, so the replication as reported tests
% a claim the paper no longer makes as its main one. This re-runs it on the
% scale now in use.
%
% Two things change, and both favour the replication having something to say.
% The binary contrast rested on the 18 hippocampal contacts that cleared the
% excursion gate; the continuous outcome is defined for every hippocampal
% contact, gate or no gate. And the spatial gradient becomes testable, which it
% was not before.
%
% Coordinates in this cohort are in native scanner space, which the dataset
% documentation states is not comparable across patients. That rules out a
% between-patient distance analysis but not a within-shaft one: distance along
% a single electrode is a within-patient quantity, and the discovery cohort's
% gradient result is itself a within-shaft effect. So the comparison is
% like-for-like on exactly the analysis that carried the discovery finding.
%
% The hippocampal reference point is each participant's own hippocampal
% centroid, computed in their own space, for the same reason.

arguments
    cfg (1,1) struct
    options.measuresFile (1,1) string = ...
        fullfile(cfg.tableDir, "replication_contact_measures_gaze_regressed_ruleA.csv")
end

if ~isfile(options.measuresFile)
    error('phg:MissingReplicationMeasures', ...
        'Replication measures not found: %s', options.measuresFile);
end
measures = readtable(options.measuresFile, 'TextType', 'string');

if ~ismember('FitHeight', measures.Properties.VariableNames)
    error('phg:ReplicationLacksContinuousOutcome', ...
        ['%s predates the continuous outcome. Re-run ' ...
         'run_replication_stage so the Gaussian amplitude is exported.'], ...
        options.measuresFile);
end

warningState = warning('off', 'all');

region = lower(strtrim(measures.NMM));
contactClass = repmat("Extrahippocampal", height(measures), 1);
contactClass(region == "hippocampus") = "Hippocampal";

scale = median(abs(measures.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end

frame = table(categorical(measures.PtID), ...
    categorical(measures.PtID + "/" + measures.Shaft), ...
    asinh(measures.FitHeight ./ scale), ...
    categorical(contactClass, ["Extrahippocampal", "Hippocampal"]), ...
    measures.FitR2, 'VariableNames', ...
    {'PtID', 'Lead', 'Coupling', 'ContactClass', 'FitR2'});
if ismember('XYZMNI_1', measures.Properties.VariableNames)
    frame.XYZ = [measures.XYZMNI_1, measures.XYZMNI_2, measures.XYZMNI_3];
elseif ismember('XYZMNI', measures.Properties.VariableNames)
    frame.XYZ = measures.XYZMNI;
else
    frame.XYZ = nan(height(measures), 3);
end
frame = frame(isfinite(frame.Coupling), :);

nHippocampal = sum(frame.ContactClass == "Hippocampal");
fprintf(['[PHG] Replication, continuous outcome: %d contacts, %d ', ...
    'hippocampal, %d participants\n'], height(frame), nHippocampal, ...
    numel(categories(removecats(frame.PtID))));

%% -------------------------------------------- Primary contrast
primary = table();
if nHippocampal >= 5 && numel(unique(frame.ContactClass)) > 1
    model = fitlme(frame, ...
        'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
        'FitMethod', 'REML');
    [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
    index = find(startsWith(string(stats.Name), "ContactClass"), 1);
    primary = table(height(frame), nHippocampal, stats.Estimate(index), ...
        stats.Lower(index), stats.Upper(index), stats.tStat(index), ...
        stats.DF(index), stats.pValue(index), ...
        stats.Estimate(index) / std(frame.Coupling), ...
        'VariableNames', {'n_contacts', 'n_hippocampal', 'estimate', ...
        'ci95_low', 'ci95_high', 't_statistic', 'satterthwaite_df', ...
        'p_value', 'standardised'});
    phg.writeTableAtomic(primary, ...
        fullfile(cfg.tableDir, 'replication_continuous_primary.csv'));
    fprintf(['[PHG]   contrast: beta = %+.3f [%+.3f, %+.3f], P = %.3g, ', ...
        '%.2f SD\n'], primary.estimate, primary.ci95_low, ...
        primary.ci95_high, primary.p_value, primary.standardised);
else
    fprintf('[PHG]   contrast not estimable: %d hippocampal contacts\n', ...
        nHippocampal);
end

%% ---------------------------------------- Within-shaft gradient
gradient = table();
hasCoordinates = all(isfinite(frame.XYZ), 2);
if sum(hasCoordinates) > 50
    spatial = frame(hasCoordinates, :);
    % Hippocampal centroid per participant, in that participant's own space.
    patients = categories(removecats(spatial.PtID));
    distance = nan(height(spatial), 1);
    for k = 1:numel(patients)
        rows = spatial.PtID == patients{k};
        hippRows = rows & spatial.ContactClass == "Hippocampal";
        if ~any(hippRows)
            continue
        end
        centroid = mean(spatial.XYZ(hippRows, :), 1);
        distance(rows) = sqrt(sum((spatial.XYZ(rows, :) - centroid) .^ 2, 2));
    end
    spatial.Distance = distance;
    spatial = spatial(isfinite(spatial.Distance), :);

    if height(spatial) > 50
        [shaftId, ~] = findgroups(string(spatial.Lead));
        shaftMean = splitapply(@mean, spatial.Distance, shaftId);
        spatial.WithinShaft = spatial.Distance - shaftMean(shaftId);
        spans = splitapply(@range, spatial.Distance, shaftId);
        spatial = spatial(spans(shaftId) > cfg.gradient.minimumShaftSpanMm, :);
        if height(spatial) > 50 && std(spatial.WithinShaft) > 0
            spatial.WithinShaftZ = spatial.WithinShaft ./ std(spatial.WithinShaft);
            model = fitlme(spatial, ...
                'Coupling ~ WithinShaftZ + (1|PtID) + (1|PtID:Lead)', ...
                'FitMethod', 'REML');
            [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
            index = find(string(stats.Name) == "WithinShaftZ", 1);
            gradient = table(height(spatial), stats.Estimate(index), ...
                stats.Lower(index), stats.Upper(index), ...
                stats.pValue(index), ...
                stats.Estimate(index) / std(spatial.Coupling), ...
                'VariableNames', {'n_contacts', 'estimate', 'ci95_low', ...
                'ci95_high', 'p_value', 'standardised'});
            phg.writeTableAtomic(gradient, ...
                fullfile(cfg.tableDir, 'replication_continuous_gradient.csv'));
            fprintf(['[PHG]   within-shaft gradient: beta = %+.3f ', ...
                '[%+.3f, %+.3f], P = %.3g, %.2f SD, n = %d\n'], ...
                gradient.estimate, gradient.ci95_low, gradient.ci95_high, ...
                gradient.p_value, gradient.standardised, gradient.n_contacts);
        end
    end
end

warning(warningState);
outputs = struct('frame', frame, 'primary', primary, 'gradient', gradient);
end
