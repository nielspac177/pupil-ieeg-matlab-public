function colors = regionColor(region)
%REGIONCOLOR Colorblind-conscious colors shared across figures.

region = string(region(:));
colors = repmat([0.62 0.62 0.62], numel(region), 1);
mapping = {
    "MTG middle temporal gyrus", [0.12 0.31 0.50]
    "Hippocampus", [0.77 0.27 0.21]
    "PrG precentral gyrus", [0.06 0.44 0.38]
    "PoG postcentral gyrus", [0.31 0.64 0.82]
    "ACgG anterior cingulate gyrus", [0.90 0.62 0.00]
    "TTG transverse temporal gyrus", [0.51 0.32 0.90]
    "MFG middle frontal gyrus", [0.42 0.42 0.42]
    "PHG parahippocampal gyrus", [0.66 0.35 0.33]
    "Thalamus Proper", [0.40 0.26 0.15]
    "Amygdala", [0.70 0.43 0.00]
    "STG superior temporal gyrus", [0.33 0.42 0.48]
    "SFG superior frontal gyrus", [0.58 0.58 0.58]
    "ITG inferior temporal gyrus", [0.48 0.55 0.58]
    "MFC medial frontal cortex", [0.55 0.55 0.55]
    "PIns posterior insula", [0.65 0.50 0.00]
    };
for k = 1:size(mapping, 1)
    colors(region == mapping{k,1}, :) = repmat(mapping{k,2}, ...
        sum(region == mapping{k,1}), 1);
end
end
