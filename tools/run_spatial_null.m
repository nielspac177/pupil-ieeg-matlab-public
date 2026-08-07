%% Spatial-null test for the receptor alignment.
repoRoot = '/Users/nielspacheco/Desktop/Research/Rolston lab/Pupils/pupil-ieeg-matlab';
addpath(fullfile(repoRoot,'src')); addpath(fullfile(repoRoot,'config'));
cfg = default_config(repoRoot);
lc = fullfile(repoRoot,'config','local_config.m'); if isfile(lc), run(lc); end
[tables,~] = phg.loadDerivedTables(cfg);
T = tables.channel;

petDir = fullfile(cfg.workDir,'cache','pet');
maps = {'NET', fullfile(petDir,'NAT_MRB_hc77_ding.nii'); ...
        'VAChT', fullfile(petDir,'VAChT_feobv_hc18_aghourian_sum.nii')};

mni = T.XYZMNI;
sampled = struct;
for m = 1:size(maps,1)
    info = niftiinfo(maps{m,2}); vol = double(niftiread(maps{m,2}));
    A = info.Transform.T';
    ijk = round(A \ [mni'; ones(1,height(T))])'; ijk = ijk(:,1:3)+1;
    ok = all(ijk>=1,2) & ijk(:,1)<=info.ImageSize(1) & ...
         ijk(:,2)<=info.ImageSize(2) & ijk(:,3)<=info.ImageSize(3);
    v = nan(height(T),1);
    v(ok) = vol(sub2ind(info.ImageSize, ijk(ok,1), ijk(ok,2), ijk(ok,3)));
    sampled.(maps{m,1}) = (v - mean(v,'omitnan'))./std(v,'omitnan');
end

region = phg.cleanRegionLabels(T.NMM);
N_SURROGATES = 500;
rows = {};

for m = 1:size(maps,1)
    marker = maps{m,1};
    for scope = ["All excursion contacts", "Neocortex only"]
        keep = T.RespAreaNet ~= 0 & isfinite(sampled.(marker));
        if scope == "Neocortex only"
            keep = keep & ~ismember(region, ["Hippocampus","Amygdala"]);
        end
        frame = table;
        frame.PtID = removecats(categorical(string(T.PtID(keep))));
        frame.Lead = categorical(phg.parseLeadLabel(string(T.Label(keep))));
        frame.IsDilation = double(T.RespAreaNet(keep) > 0);
        frame.Value = sampled.(marker)(keep);
        coords = mni(keep,:);

        formula = 'IsDilation ~ Value + (1|PtID) + (1|PtID:Lead)';
        model = fitglme(frame, formula, 'Distribution','Binomial', ...
            'Link','logit','FitMethod','Laplace');
        observed = model.Coefficients.Estimate(strcmp(model.Coefficients.Name,'Value'));
        parametricP = model.Coefficients.pValue(strcmp(model.Coefficients.Name,'Value'));

        S = phg.spatialNullSurrogates(frame.Value, coords, N_SURROGATES);
        nullBeta = nan(N_SURROGATES,1);
        surrogateFrame = frame;
        for k = 1:N_SURROGATES
            surrogateFrame.Value = S(:,k);
            try
                nullModel = fitglme(surrogateFrame, formula, ...
                    'Distribution','Binomial','Link','logit','FitMethod','Laplace');
                nullBeta(k) = nullModel.Coefficients.Estimate( ...
                    strcmp(nullModel.Coefficients.Name,'Value'));
            catch
            end
        end
        valid = isfinite(nullBeta);
        spatialP = (1 + sum(abs(nullBeta(valid)) >= abs(observed))) / (1 + sum(valid));

        fprintf(['[PHG] %-6s %-22s beta %+.3f (OR %.2f)  parametric P = %.3g  ' ...
            'spatial-null P = %.3f  (%d surrogates)\n'], marker, scope, ...
            observed, exp(observed), parametricP, spatialP, sum(valid));

        rows{end+1,1} = table(string(marker), string(scope), height(frame), ...
            exp(observed), parametricP, spatialP, sum(valid), ...
            std(nullBeta(valid)), ...
            'VariableNames', {'marker','sample','n_contacts','odds_ratio', ...
            'parametric_p_value','spatial_null_p_value','n_surrogates', ...
            'null_beta_sd'}); %#ok<SAGROW>
    end
end

result = vertcat(rows{:});
phg.writeTableAtomic(result, fullfile(cfg.tableDir,'receptor_spatial_null.csv'));
disp(result);
fprintf('[PHG] Spatial-null test complete.\n');
