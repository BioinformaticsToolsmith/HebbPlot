% Author: Alfredo Velasco
% This script will find the average of the matrix and create its hebbPlot
% (except it's not a HEBBPlot because no Hebbian learning is used)

A = dlmread('/Volumes/Bioinf/ExpPE003/results/matrix.txt');
markList = readMarks('/Volumes/Bioinf/ExpPE003/results/marks.txt', 0);

m = mean(A);
[s, markIndex] = sortASampleHierarchical(m, 27, 41, 5);
% F = plotASample( s, myR, myC, marks(markIndex), colorMapIn );
F = plotASample(s, 27, 41, markList(markIndex), jet, [min(m) max(m)])
saveAsPdf(F, '/Users/alfredo/Writing/HebbPlot/LaTexBMC/Figures/E003Individual/averageMatrix.pdf');