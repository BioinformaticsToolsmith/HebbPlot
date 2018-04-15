% Author: Alfredo Velasco
%  This is the hierarchical clustering experiment.
%  What this does is take in a matrix file and convert it from its
%  column-wise format to a row-wise property. It will then print out the
%  matrix into a heatmap, just like HebbPlot.

%Read a matrix
% Reshape each row (filled column-wise) --> matrix --> row (filled row-wise)
% --> end up with another matrix with the same dimension
% Run hierarchical clustering on this matrix

matrixFile = '/Volumes/Bioinf/ExpPE003/results/matrix.txt';
markFile = '/Volumes/Bioinf/ExpPE003/results/marks.txt';
pdfFile = '/Users/alfredo/Writing/HebbPlot/LaTexBMC/Figures/E003Individual/hierarchicalH1P300.pdf';
% file = '/Users/alfredo/Projects/Histone/Perl/HebbPlotOriginal/testMatrix.txt';
W = dlmread(matrixFile);
% W = W(1:10, :)
% W = round(rand(size(W))) * 2 - 1 ;
% N will be the matrix file as a row-wise representation of the overlaps
N = zeros(size(W));
[r, c] = size(W);
for i = 1:r
    cWise = W(i, :);
    cWiseMatrix = reshape(cWise, c / 41, 41);
    cWiseMatrix = cWiseMatrix';
    N(i, :) = transpose(cWiseMatrix(:));
end

clear W;
% nCopy = N;
% for i = 4:r-4
%     for j = 4:c-4
%         N(i, j) = mean(nCopy(i,j-3:j+3));
%     end
% end
% clear nCopy;
% N = N + rand(size(N));
% 1

% This part does the clustering of the matrix N
cityBlock = 'cityblock';
Z = linkage(N, 'weighted', cityBlock);
D = pdist(N, cityBlock);
enhancerIndex = optimalleaforder(Z, D);

% Generate the figure
h = figure;
imagesc(N(enhancerIndex, :));colormap(jet);

colorbar;
markLabels = readMarks(markFile , 0)
center = size(markLabels, 2)
interval = round(c / (center))
%interval:interval:c-interval

t = 1:interval:c;
t = t + round(interval / 2)
set(gca, 'XTick', t , 'XTickLabel', markLabels);
set(gca, 'FontName', 'Times');
xtickangle(60);
set(gca, 'YTick', []);

xlabel('Marks', 'interpreter', 'latex', 'FontName', 'Times');
ylabel('Enhancers', 'interpreter', 'latex', 'FontName', 'Times');

% Set the font size of all text in the figure
set(gca, 'FontSize', 8);

%Need to save as PDF
saveAsPdf(h, pdfFile);