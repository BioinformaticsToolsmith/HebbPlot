function processTwoExp( expDir , othExpDir, vLine, ex)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%PROCESSTWOEXP sorts Hebb Plots of the second experiment using the order of
%the first

l = dir(expDir);

for i=3:size(l,1)
    % Matrix file
    [w, markIndex, markList] = helper(l(i), expDir, vLine, ex);
    [w1, ignore1, ignore2] = helper(l(i), othExpDir, vLine, ex);
    s = reshape(w1, size(w1, 2)/ vLine, vLine);
    F = plotASample(s(markIndex, :), size(w1, 2)/ vLine, vLine, ...
        markList(markIndex));
    pdfFile = [othExpDir '/' l(i).name '/HebbPlot' l(i).name '.pdf'];
    saveAsPdf(F, pdfFile);
    close all;
end
end


function [w, markIndex, markList] = helper(l, expDir, vLine, ex)
% Matrix file
matFile = [expDir '/' l.name '/results/matrix.txt'];

disp(matFile);

if(exist(matFile, 'file') == 2)
    disp(['Processing ' expDir '/' l.name '...']);
    % Read marks
    markFile = [expDir '/' l.name '/results/marks.txt'];
    id = fopen(markFile);
    markList = textscan(id, '%s', 'Delimiter' , ',');
    markList = reshape(markList{1}, 1, size(markList{1},1));
    markList = strrep(markList, '''', '');
    fclose(id);
    
    % Pdf file
    pdfFile = [expDir '/' l.name '/results/hebbPlotRearranged.pdf'];
    histFile = [expDir '/' l.name '/results/HistRearranged.pdf'];
    
    % This calculates the learning rate
    matrix = dlmread (matFile);
    [numRows,~] = size(matrix);
    learningRate = ((10000 / numRows) * 0.001);
    
    
    % Make plot
    [w, markIndex] = hebbPlot(matFile, vLine, learningRate, markList, ex, ...
        pdfFile, histFile);
    wFile = [expDir '/' l.name '/results/w.mat'];
    save(wFile, 'w');
    %close all;
else
    disp(matFile);
end
end