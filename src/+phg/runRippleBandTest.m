function outputs = runRippleBandTest(wavelet, channel, cfg)
%RUNRIPPLEBANDTEST Is hippocampal pupil coupling band-limited near the ripple range?
%
%   The mechanistic account offered for the polarity reversal is that
%   hippocampal high-gamma peaks are dominated by sharp-wave ripples, which
%   occur in quiescent low-arousal states. Testing that properly needs the
%   continuous mesial-temporal LFP, which is not in the supplied archive: a
%   ripple is a discrete 80-150 Hz event and cannot be recovered from an
%   averaged cross-correlation spectrum.
%
%   One weaker prediction *can* be tested with what exists. If hippocampal
%   coupling is carried by ripples, its frequency profile should be relatively
%   concentrated in the ripple band, whereas broadband neocortical high-gamma
%   coupling should not be. Each contact's profile is normalised by its own
%   mean over 30-256 Hz first, so the comparison is about spectral SHAPE and
%   not about which contacts couple more strongly.
%
%   This is a shape test, not a ripple detector. A positive result is
%   consistent with the ripple account and does not establish it; a null
%   result does not refute it, because averaging over peaks blurs exactly the
%   band-limited structure the test looks for. It is reported as exploratory.

arguments
    wavelet table
    channel table
    cfg (1,1) struct
end

RIPPLE_BAND = [80 150];
REFERENCE_BAND = [30 256];

wavelet.Region = phg.cleanRegionLabels(wavelet.NMM);
wavelet.ContactKey = string(wavelet.PtID) + "|" + string(wavelet.Label);

% Only contacts that showed a broadband excursion enter the test: on a contact
% with no coupling at all the normalised profile is noise divided by noise.
channel.ContactKey = string(channel.PtID) + "|" + string(channel.Label);
excursionKeys = channel.ContactKey(channel.RespAreaNet ~= 0);
wavelet = wavelet(ismember(wavelet.ContactKey, excursionKeys), :);

inReference = wavelet.freq >= REFERENCE_BAND(1) & wavelet.freq <= REFERENCE_BAND(2);
wavelet = wavelet(inReference, :);

[group, contactKey] = findgroups(wavelet.ContactKey);
nContacts = numel(contactKey);
index = nan(nContacts, 1);
for k = 1:nContacts
    rows = wavelet(group == k, :);
    amplitude = abs(rows.RespAreaAbs);
    reference = mean(amplitude, 'omitnan');
    if ~isfinite(reference) || reference <= 0
        continue
    end
    normalised = amplitude ./ reference;
    inRipple = rows.freq >= RIPPLE_BAND(1) & rows.freq <= RIPPLE_BAND(2);
    if ~any(inRipple) || all(inRipple)
        continue
    end
    % log2 ratio of in-band to out-of-band normalised amplitude: zero means the
    % ripple band is no more represented than the rest of the profile.
    index(k) = log2(mean(normalised(inRipple), 'omitnan') ./ ...
        mean(normalised(~inRipple), 'omitnan'));
end

% Rebuild the per-contact design from the channel table so that region, patient
% and shaft come from one authority.
[~, position] = ismember(contactKey, channel.ContactKey);
keep = position > 0 & isfinite(index);
contactTable = table(contactKey(keep), index(keep), ...
    categorical(string(channel.PtID(position(keep)))), ...
    categorical(phg.parseLeadLabel(string(channel.Label(position(keep))))), ...
    phg.cleanRegionLabels(channel.NMM(position(keep))), ...
    'VariableNames', {'ContactKey', 'RippleIndex', 'PtID', 'Lead', 'Region'});
contactTable.ContactClass = categorical(repmat("Extrahippocampal", ...
    height(contactTable), 1));
contactTable.ContactClass(contactTable.Region == "Hippocampus") = "Hippocampal";
contactTable.ContactClass = reordercats(contactTable.ContactClass, ...
    {'Extrahippocampal', 'Hippocampal'});

model = fitlme(contactTable, ...
    'RippleIndex ~ ContactClass + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
[~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');

term = string(stats.Name);
resultTable = table(term, stats.Estimate, stats.SE, stats.tStat, stats.DF, ...
    stats.pValue, stats.Lower, stats.Upper, ...
    'VariableNames', {'term', 'estimate_log2', 'standard_error', ...
    't_statistic', 'satterthwaite_df', 'p_value', 'ci95_low', 'ci95_high'});
phg.writeTableAtomic(resultTable, ...
    fullfile(cfg.tableDir, 'ripple_band_selectivity.csv'));

isHippocampal = contactTable.ContactClass == "Hippocampal";
summaryTable = table( ...
    ["Hippocampal"; "Extrahippocampal"], ...
    [sum(isHippocampal); sum(~isHippocampal)], ...
    [numel(unique(contactTable.PtID(isHippocampal))); ...
     numel(unique(contactTable.PtID(~isHippocampal)))], ...
    [median(contactTable.RippleIndex(isHippocampal), 'omitnan'); ...
     median(contactTable.RippleIndex(~isHippocampal), 'omitnan')], ...
    'VariableNames', {'contact_class', 'n_contacts', 'n_patients', ...
    'median_log2_ripple_index'});
phg.writeTableAtomic(summaryTable, ...
    fullfile(cfg.tableDir, 'ripple_band_summary.csv'));

outputs = struct('resultTable', resultTable, 'summaryTable', summaryTable, ...
    'contactTable', contactTable);

effect = resultTable(startsWith(resultTable.term, "ContactClass"), :);
if ~isempty(effect)
    fprintf(['[PHG] Ripple-band shape test (exploratory): hippocampal minus ' ...
        'extrahippocampal log2 index = %+.3f (P = %.3g) on %d contacts.\n'], ...
        effect.estimate_log2(1), effect.p_value(1), height(contactTable));
end
end
