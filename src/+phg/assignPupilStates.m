function [events, exposure] = assignPupilStates(events, pupilTime, pupilSize, nStates)
%ASSIGNPUPILSTATES Assign accepted ripple events to pupil-size quantiles.

arguments
    events table
    pupilTime double
    pupilSize double
    nStates (1,1) double {mustBeInteger,mustBePositive} = 6
end
pupilTime = pupilTime(:);
pupilSize = pupilSize(:);
if numel(pupilTime) ~= numel(pupilSize)
    error('phg:PupilLengthMismatch', ...
        'pupilTime and pupilSize must have equal length.');
end
valid = isfinite(pupilTime) & isfinite(pupilSize);
if sum(valid) < nStates * 10
    error('phg:InsufficientPupilData', ...
        'Not enough finite pupil samples for %d states.', nStates);
end

edges = quantile(pupilSize(valid), linspace(0, 1, nStates + 1));
edges(1) = -Inf;
edges(end) = Inf;
for k = 2:numel(edges)
    if edges(k) <= edges(k-1)
        edges(k) = edges(k-1) + eps(edges(k-1));
    end
end
sampleState = discretize(pupilSize, edges);

centerTime = (events.onset_seconds + events.offset_seconds) ./ 2;
eventPupil = interp1(pupilTime(valid), pupilSize(valid), centerTime, ...
    'linear', NaN);
events.pupil_size = eventPupil;
events.pupil_state = discretize(eventPupil, edges);

sampleInterval = median(diff(pupilTime(valid)), 'omitnan');
state = (1:nStates)';
durationSeconds = zeros(nStates,1);
for k = 1:nStates
    durationSeconds(k) = sum(sampleState == k) * sampleInterval;
end
exposure = table(state, edges(1:end-1)', edges(2:end)', durationSeconds, ...
    'VariableNames', {'pupil_state', 'lower_edge', 'upper_edge', ...
    'duration_seconds'});
end
