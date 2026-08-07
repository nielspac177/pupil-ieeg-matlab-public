function qc = exportPublicationFigure(figureHandle, cfg, stem, options)
%EXPORTPUBLICATIONFIGURE Export manuscript formats and boundary-QC report.

arguments
    figureHandle (1,1) matlab.ui.Figure
    cfg (1,1) struct
    stem (1,1) string
    options.pdfContentType (1,1) string = "vector"
    options.saveEditable (1,1) logical = true
end

drawnow;
base = fullfile(cfg.figureDir, stem);
pngPath = base + ".png";
originalUnits = figureHandle.Units;
figureHandle.Units = 'inches';
figureSize = figureHandle.Position(3:4);
figureHandle.PaperUnits = 'inches';
figureHandle.PaperPosition = [0 0 figureSize];
figureHandle.PaperSize = figureSize;
figureHandle.PaperPositionMode = 'manual';
figureHandle.InvertHardcopy = 'off';
print(figureHandle, char(pngPath), '-dpng', '-r600');
print(figureHandle, char(base + ".tiff"), '-dtiff', '-r600');
if options.pdfContentType == "image"
    print(figureHandle, char(base + ".pdf"), '-dpdf', '-image', '-r600');
else
    print(figureHandle, char(base + ".pdf"), '-dpdf', '-vector', '-r600');
end
figureHandle.Units = originalUnits;
if options.saveEditable
    savefig(figureHandle, base + ".fig");
end
qc = phg.auditFigureLayout(figureHandle, pngPath, stem);
phg.writeTableAtomic(qc, fullfile(cfg.tableDir, ...
    "figure_qc_" + stem + ".csv"));
if ~qc.boundary_clear
    warning('phg:FigureBoundary', ...
        'Figure %s has content within four pixels of an export boundary.', stem);
end
if ~qc.margins_tight
    warning('phg:FigureDeadMargin', ...
        'Figure %s has a dead margin of %.2f cm at its published size.', ...
        stem, qc.largest_margin_cm);
end
if qc.text_overlap_count > 0
    warning('phg:FigureTextOverlap', ...
        'Figure %s has %d overlapping text boxes (worst %.0f%%: %s).', ...
        stem, qc.text_overlap_count, 100*qc.worst_overlap_fraction, ...
        qc.worst_overlap_pair);
end
end
