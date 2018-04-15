function saveAsPdf(h, outFile)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
% Saves the figure handler h as a pdf.
% The pdf will be the outFile.
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,outFile,'-dpdf','-r0')
