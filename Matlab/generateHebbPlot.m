

function [w, markIndex] = generateHebbPlot(inOutDir, vLines, eLines, colorMapIn)
%
% Author: Hani Zakaria Girgis, PhD
% Author: Alfredo Velasco II
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
% GENERATEHEBBPLOT this is the driver to generate a hebbPlot
% inOutDir: is a directory containing the matrix file, the marks file, and
% the regions file. Upon successful completion, this directory will
% include the HebbPlot, the quality histogram, and the weight vector. 
% vLines: is the number of uniform samples per region
% eLines: is the number of lines to be excluded on both sides, must be > 0

matrixFile = [inOutDir '/matrix.txt'];
marksFile = [inOutDir '/marks.txt'];
hebbFile = [inOutDir '/hebbPlot.pdf'];
histFile = [inOutDir '/qualityHist.pdf'];
wFile = [inOutDir '/w.mat'];

id = fopen(marksFile);
markList = textscan(id, '%s', 'Delimiter' , ',');
markList = reshape(markList{1}, 1, size(markList{1},1));
fclose(id);
readFile = dlmread (matrixFile);
[numRows,~] = size(readFile);
learningRate = ((10000 / numRows) * 0.001);

% hebbPlot(matrixFile, vLines, 0.001, markList, eLines, hebbFile, histFile);
[w, markIndex] = hebbPlot(matrixFile, vLines, learningRate, markList, eLines, ...
    colorMapIn ,hebbFile, histFile);
save(wFile, 'w');
end