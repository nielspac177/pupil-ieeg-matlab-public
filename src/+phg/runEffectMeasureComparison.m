function outputs = runEffectMeasureComparison(T, cfg)
%RUNEFFECTMEASURECOMPARISON Report the primary contrast on three scales.
%
% The primary model is a binomial GLMM, so its natural parameter is an odds
% ratio. That is a defensible model parameter and a poor effect size here,
% because the outcome is common: 69% of excursion contacts are dilation-linked.
% The odds ratio and the risk ratio coincide only for rare outcomes, and
% diverge further as prevalence rises. At this prevalence the divergence is
% fivefold, so "odds ratio 0.065" and "a fifteen-fold reduction" describe a
% threefold reduction in probability.
%
% This is not a subtlety about wording. A reader who takes an odds ratio for a
% risk ratio overstates the effect by a factor of five.
%
% The deeper issue is that this field does not use odds ratios at all. The two
% closest precedents both model a *continuous signed* coupling value rather
% than binarising its sign: Alasfour et al. (2021) compute a Spearman
% correlation per contact, Fisher z-transform it and fit a linear mixed model
% with region as a predictor, stating explicitly that they do not binarise;
% Pfeffer et al. (2022) use Spearman correlations, mutual information and
% cross-correlation lags with FDR-corrected t-tests. Neither reports an odds
% ratio. Binarising the sign is therefore both non-standard here and lossy,
% since it discards the magnitude the field actually models.
%
% See docs/effect_measure_review.md for what follows from that.
%
% Three measures are therefore exported together:
%   odds ratio       the model parameter, from the logistic GLMM
%   risk ratio       from a log-binomial GLMM, cross-checked with a
%                    modified-Poisson fit because log-binomial models
%                    frequently fail to converge
%   risk difference  absolute, in percentage points, which needs no scale
%
% The manuscript should lead with the absolute proportions and the risk
% difference, and name the odds ratio as a model parameter rather than as an
% effect size.

arguments
    T table
    cfg (1,1) struct
end

region = phg.cleanRegionLabels(T.NMM);
isHippocampus = region == "Hippocampus";
contactClass = repmat("Extrahippocampal", height(T), 1);
contactClass(isHippocampus) = "Hippocampal";

frame = table(categorical(string(T.PtID)), ...
    categorical(string(T.PtID) + "/" + phg.parseLeadLabel(T.Label)), ...
    double(T.RespAreaNet > 0), ...
    categorical(contactClass, ["Extrahippocampal", "Hippocampal"]), ...
    'VariableNames', {'PtID', 'Lead', 'IsDilation', 'ContactClass'});
frame = frame(T.RespAreaNet ~= 0, :);

hippocampal = frame.ContactClass == "Hippocampal";
riskHippocampal = mean(frame.IsDilation(hippocampal));
riskOther = mean(frame.IsDilation(~hippocampal));
prevalence = mean(frame.IsDilation);

formula = 'IsDilation ~ ContactClass + (1|PtID) + (1|PtID:Lead)';
oddsModel = fitglme(frame, formula, 'Distribution', 'Binomial', ...
    'Link', 'logit', 'FitMethod', 'Laplace');
[oddsRatio, oddsLow, oddsHigh, oddsP] = localTerm(oddsModel, "ContactClass");

% Log-binomial is the textbook estimator and the one that most often fails to
% converge, because nothing constrains its fitted probabilities below one.
% Modified Poisson is fitted regardless and reported beside it; agreement
% between the two is the check that neither is an artefact of its own
% parameterisation.
riskRatio = NaN; riskLow = NaN; riskHigh = NaN; riskP = NaN;
riskConverged = false;
warningState = warning('off', 'all');
try
    riskModel = fitglme(frame, formula, 'Distribution', 'Binomial', ...
        'Link', 'log', 'FitMethod', 'Laplace');
    [riskRatio, riskLow, riskHigh, riskP] = localTerm(riskModel, "ContactClass");
    riskConverged = isfinite(riskRatio);
catch
    riskConverged = false;
end
poissonModel = fitglme(frame, formula, 'Distribution', 'Poisson', ...
    'Link', 'log', 'FitMethod', 'Laplace');
[poissonRatio, poissonLow, poissonHigh, poissonP] = ...
    localTerm(poissonModel, "ContactClass");
warning(warningState);

measure = ["odds_ratio_logistic_glmm"
           "risk_ratio_log_binomial_glmm"
           "risk_ratio_modified_poisson"
           "risk_difference_percentage_points"];
estimate = [oddsRatio; riskRatio; poissonRatio; ...
    100 * (riskHippocampal - riskOther)];
ciLow = [oddsLow; riskLow; poissonLow; NaN];
ciHigh = [oddsHigh; riskHigh; poissonHigh; NaN];
pValue = [oddsP; riskP; poissonP; NaN];
converged = [true; riskConverged; true; true];

comparison = table(measure, estimate, ciLow, ciHigh, pValue, converged, ...
    'VariableNames', {'measure', 'estimate', 'ci95_low', 'ci95_high', ...
    'p_value', 'converged'});
phg.writeTableAtomic(comparison, ...
    fullfile(cfg.tableDir, 'primary_effect_measures.csv'));

context = table(prevalence, riskHippocampal, riskOther, ...
    (1 / oddsRatio) / (1 / poissonRatio), ...
    'VariableNames', {'outcome_prevalence', 'risk_hippocampal', ...
    'risk_extrahippocampal', 'odds_to_risk_ratio_inflation'});
phg.writeTableAtomic(context, ...
    fullfile(cfg.tableDir, 'primary_effect_measure_context.csv'));

fprintf(['[PHG] Effect measures: OR %.3f, RR %.3f (Poisson %.3f), ', ...
    'RD %.1f points; outcome prevalence %.0f%%, OR inflated %.1fx\n'], ...
    oddsRatio, riskRatio, poissonRatio, ...
    100 * (riskHippocampal - riskOther), 100 * prevalence, ...
    (1 / oddsRatio) / (1 / poissonRatio));

outputs = struct('comparison', comparison, 'context', context);
end

% -------------------------------------------------------------------------
function [ratio, low, high, p] = localTerm(model, prefix)
stats = model.Coefficients;
index = find(startsWith(string(stats.Name), prefix), 1);
if isempty(index)
    [ratio, low, high, p] = deal(NaN);
    return
end
ratio = exp(stats.Estimate(index));
low = exp(stats.Lower(index));
high = exp(stats.Upper(index));
p = stats.pValue(index);
end
