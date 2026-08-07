function style = figureStyle()
%FIGURESTYLE Shared journal-scale typography and color definitions.

style = struct;
style.fontName = 'Arial';
style.tickFontSize = 7;
style.labelFontSize = 8;
style.titleFontSize = 8;
style.panelFontSize = 10;
style.lineWidth = 1.15;
style.axisLineWidth = 0.75;
style.gray = [0.58 0.60 0.63];
style.lightGray = [0.78 0.80 0.82];
style.dilation = [0.90 0.35 0.24];
style.constriction = [0.15 0.29 0.34];
style.sensorimotor = [0.10 0.49 0.42];
style.task = [0.10 0.49 0.42; 0.90 0.62 0.00; 0.12 0.31 0.50];
end
