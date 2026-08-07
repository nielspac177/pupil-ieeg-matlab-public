function label = shortRegionLabel(region)
%SHORTREGIONLABEL Compact anatomical labels for publication axes.

region = string(region);
label = region;
mapping = {
    "MTG middle temporal gyrus", "MTG"
    "Hippocampus", "Hipp"
    "PrG precentral gyrus", "PrG"
    "PoG postcentral gyrus", "PoG"
    "ACgG anterior cingulate gyrus", "ACgG"
    "TTG transverse temporal gyrus", "TTG"
    "MFG middle frontal gyrus", "MFG"
    "PHG parahippocampal gyrus", "PHG"
    "Thalamus Proper", "Thal"
    "Amygdala", "Amy"
    "STG superior temporal gyrus", "STG"
    "SFG superior frontal gyrus", "SFG"
    "ITG inferior temporal gyrus", "ITG"
    "MFC medial frontal cortex", "MFC"
    "PIns posterior insula", "pIns"
    "OrIFG orbital part of the inferior frontal gyrus", "OrIFG"
    "TrIFG triangular part of the inferior frontal gyrus", "TrIFG"
    "POrG posterior orbital gyrus", "POrG"
    "MOrG medial orbital gyrus", "MOrG"
    "AIns anterior insula", "aIns"
    "PP planum polare", "PP"
    "CO central operculum", "CO"
    "PT planum temporale", "PT"
    "TMP temporal pole", "TMP"
    "GRe gyrus rectus", "GRe"
    "LOrG lateral orbital gyrus", "LOrG"
    };
for k = 1:size(mapping, 1)
    label(region == mapping{k,1}) = mapping{k,2};
end
long = strlength(label) > 12;
label(long) = extractBefore(label(long), 13);
end
