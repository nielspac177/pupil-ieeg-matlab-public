function decision = replicationVerdict(oddsRatio, ciLow, ciHigh, discoveryInterval)
%REPLICATIONVERDICT Apply the prespecified H1 decision rule.
%   Fixed in docs/replication_plan_ebrains.md section 2 before the data
%   arrived, and lifted out of the analysis driver so that it can be unit
%   tested against the wording rather than trusted.
%
%   "replication"                            OR < 1, 95% interval excludes 1
%   "directionally_consistent_underpowered"  OR < 1, interval includes 1, and
%                                            OR lies inside the discovery
%                                            interval [0.021, 0.207]
%   "inconclusive"                           OR < 1, interval includes 1, OR
%                                            outside the discovery interval
%   "failure"                                OR >= 1, or an interval excluding
%                                            1 in the opposite direction
%   "not_estimable"                          the model could not be fitted
%
%   The parent protocol is explicit that the middle category is *not*
%   replication, and this function refuses to blur them.

arguments
    oddsRatio (1,1) double
    ciLow (1,1) double
    ciHigh (1,1) double
    discoveryInterval (1,2) double
end

if ~isfinite(oddsRatio) || ~isfinite(ciLow) || ~isfinite(ciHigh)
    decision = "not_estimable";
elseif oddsRatio >= 1
    decision = "failure";
elseif ciLow > 1
    decision = "failure";
elseif ciHigh < 1
    decision = "replication";
elseif oddsRatio >= discoveryInterval(1) && oddsRatio <= discoveryInterval(2)
    decision = "directionally_consistent_underpowered";
else
    decision = "inconclusive";
end
end
