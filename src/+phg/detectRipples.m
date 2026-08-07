function events = detectRipples(signalUv, samplingFrequency, options)
%DETECTRIPPLES Detect candidate ripple-band events in continuous iEEG.
% SIGNALUV is samples-by-channels in microvolts. The detector returns all
% candidates and explicit rejection reasons so thresholds can be audited.

arguments
    signalUv double
    samplingFrequency (1,1) double {mustBePositive}
    options.bandHz (1,2) double = [80 120]
    options.filterOrder (1,1) double {mustBeInteger,mustBePositive} = 4
    options.rmsWindowSeconds (1,1) double {mustBePositive} = 0.020
    options.minimumZ (1,1) double = 2.5
    options.maximumZ (1,1) double = 9.0
    options.minimumDurationSeconds (1,1) double {mustBePositive} = 0.038
    options.maximumRawAmplitudeUv (1,1) double {mustBePositive} = 300
    options.spectralWindowSeconds (1,1) double {mustBePositive} = 0.250
end

if isvector(signalUv)
    signalUv = signalUv(:);
end
if samplingFrequency <= 2 * options.bandHz(2)
    error('phg:SamplingTooLow', ...
        'Sampling frequency must exceed twice the upper ripple frequency.');
end
if options.minimumZ >= options.maximumZ
    error('phg:InvalidRippleThresholds', ...
        'minimumZ must be lower than maximumZ.');
end

nSamples = height(signalUv);
nChannels = width(signalUv);
windowSamples = max(1, round(options.rmsWindowSeconds * samplingFrequency));
minimumSamples = max(1, ceil(options.minimumDurationSeconds * samplingFrequency));

[b, a] = butter(options.filterOrder, ...
    options.bandHz ./ (samplingFrequency / 2), 'bandpass');

channelColumn = zeros(0,1);
onsetSample = zeros(0,1);
offsetSample = zeros(0,1);
durationSeconds = zeros(0,1);
peakZ = zeros(0,1);
peakFrequencyHz = zeros(0,1);
rawPeakUv = zeros(0,1);
accepted = false(0,1);
reason = strings(0,1);

for channel = 1:nChannels
    raw = signalUv(:, channel);
    finite = isfinite(raw);
    if sum(finite) < max(10 * minimumSamples, round(samplingFrequency))
        continue
    end
    filled = fillmissing(raw, 'linear', 'EndValues', 'nearest');
    filtered = filtfilt(b, a, filled);
    envelope = sqrt(movmean(filtered .^ 2, windowSamples, ...
        'Endpoints', 'shrink'));
    baselineMean = mean(envelope(finite), 'omitnan');
    baselineSd = std(envelope(finite), 'omitnan');
    if ~isfinite(baselineSd) || baselineSd <= eps
        continue
    end
    zEnvelope = (envelope - baselineMean) ./ baselineSd;
    above = zEnvelope >= options.minimumZ & finite;
    boundaries = diff([false; above; false]);
    starts = find(boundaries == 1);
    stops = find(boundaries == -1) - 1;

    for candidate = 1:numel(starts)
        first = starts(candidate);
        last = stops(candidate);
        localZ = max(zEnvelope(first:last));
        localRawPeak = max(abs(filled(first:last)));
        localDuration = (last - first + 1) / samplingFrequency;
        localReason = "accepted";
        localAccepted = true;

        if (last - first + 1) < minimumSamples
            localReason = "too_short";
            localAccepted = false;
        elseif localZ > options.maximumZ
            localReason = "above_maximum_z";
            localAccepted = false;
        elseif localRawPeak > options.maximumRawAmplitudeUv
            localReason = "raw_amplitude_artifact";
            localAccepted = false;
        elseif any(~finite(first:last))
            localReason = "overlaps_missing_data";
            localAccepted = false;
        end

        [spectralOk, frequency] = spectralVerification( ...
            filled, first, last, samplingFrequency, options);
        if localAccepted && ~spectralOk
            localReason = "failed_spectral_verification";
            localAccepted = false;
        end

        channelColumn(end+1,1) = channel; %#ok<AGROW>
        onsetSample(end+1,1) = first; %#ok<AGROW>
        offsetSample(end+1,1) = last; %#ok<AGROW>
        durationSeconds(end+1,1) = localDuration; %#ok<AGROW>
        peakZ(end+1,1) = localZ; %#ok<AGROW>
        peakFrequencyHz(end+1,1) = frequency; %#ok<AGROW>
        rawPeakUv(end+1,1) = localRawPeak; %#ok<AGROW>
        accepted(end+1,1) = localAccepted; %#ok<AGROW>
        reason(end+1,1) = localReason; %#ok<AGROW>
    end
end

eventId = (1:numel(channelColumn))';
onsetSeconds = (onsetSample - 1) ./ samplingFrequency;
offsetSeconds = (offsetSample - 1) ./ samplingFrequency;
events = table(eventId, channelColumn, onsetSample, offsetSample, ...
    onsetSeconds, offsetSeconds, durationSeconds, peakZ, peakFrequencyHz, ...
    rawPeakUv, accepted, reason, 'VariableNames', {'event_id', 'channel', ...
    'onset_sample', 'offset_sample', 'onset_seconds', 'offset_seconds', ...
    'duration_seconds', 'peak_z', 'peak_frequency_hz', 'raw_peak_uv', ...
    'accepted', 'rejection_reason'});
end

function [ok, peakFrequency] = spectralVerification( ...
    signal, first, last, samplingFrequency, options)
center = round((first + last) / 2);
halfWindow = round(options.spectralWindowSeconds * samplingFrequency / 2);
index = max(1, center-halfWindow):min(numel(signal), center+halfWindow);
segment = detrend(signal(index));
if numel(segment) < round(0.1 * samplingFrequency)
    ok = false;
    peakFrequency = NaN;
    return
end

welchLength = min(numel(segment), max(64, round(0.1 * samplingFrequency)));
nfft = max(256, 2^nextpow2(numel(segment)));
[power, frequency] = pwelch(segment, hamming(welchLength), ...
    floor(welchLength/2), nfft, samplingFrequency);
band = frequency >= options.bandHz(1) & frequency <= options.bandHz(2);
if ~any(band)
    ok = false;
    peakFrequency = NaN;
    return
end

bandPower = power(band);
bandFrequency = frequency(band);
[maximumPower, maximumIndex] = max(bandPower);
peakFrequency = bandFrequency(maximumIndex);
if maximumPower <= 0 || ~isfinite(maximumPower)
    ok = false;
    return
end

[peaks, ~] = findpeaks(bandPower, ...
    'MinPeakProminence', 0.10 * maximumPower);
substantialPeaks = peaks >= 0.80 * maximumPower;
ok = sum(substantialPeaks) == 1;
end
