

function [w, markIndex] = hebbPlot(fileName, sampleNum, alpha, marks, e, ...
    colorMapIn ,outFile, histFile)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%
% e is the number of columns to exclude on each side while sorting the 
% matrix. All columns will be displayed regardless of e. It only affects
% the sorting algorithm.
%

% Read data
d = dlmread(fileName);



% Extract samples without the identifiers
tmp = d(:,1:end);

% Denoising by removing really black rows
black = -1 * ones(1,size(tmp,2));

dSimList = zeros(size(tmp,1),1);
for i=1:size(tmp, 1)
   dSimList(i) = dotSim(black, tmp(i, :));
end

% P0 has no black vectors, which are noise
p0 = tmp(dSimList < 0.8, :); 


p = ones(size(p0));

% Make triplets
r = size(p0, 1);
p0_1 = p0(randperm(r), :);
p0_2 = p0(randperm(r), :);
p0_3 = p0;
p0_3( p0_1 ~= p0) = 0;
p0_3( p0_2 ~= p0) = 0;

% Make the prototype
w = trainOutstar( p0_3, p, alpha );

% Evaluate the quality of the prototype
h = zeros(1,r);
for j=1:r
    h(1,j) = dotSim(w, p0_3(j, :));
end
good = round(100 * size(find(h > 0.5), 2)/r);
disp(['Samples with dotsim above 0.5 = ' num2str(good) '%']);
Z = figure; 
hist(h, [-0.9:0.1:0.9]);
saveAsPdf(Z, histFile);

% Make the Hebb Plot
myR = size(p0,2) / sampleNum;
myC = sampleNum;

[s, markIndex] = sortASampleHierarchical(w, myR, myC, e);
F = plotASample( s, myR, myC, marks(markIndex), colorMapIn );
saveAsPdf(F, outFile);

% This better work
for i =1:10:100
    hani = reshape(p0(i, :), myR, myC);
    hani = hani(markIndex, :);
    hani = reshape(hani, 1, myR * myC);
    G = plotASample(hani, myR, myC, marks(markIndex), colorMapIn);
    newOutFile = strrep(outFile, '.pdf', strcat(num2str(i), '.pdf'));
    saveAsPdf(G, newOutFile);
end
end