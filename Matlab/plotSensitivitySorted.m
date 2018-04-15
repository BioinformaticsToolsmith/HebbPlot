% numList = [200 400 600 800 1000 1500 2000 2500 3000 3500 4000 4500 5000]

function plotSensitivitySorted( expDir, exDir1, exDir2, vLine, ex, numList )
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%PROCESSTWOEXP sorts Hebb Plots of the second experiment using the order of
%the first
[w, markIndex, markList] = helper([ exDir1 num2str(numList(1)) exDir2 ], expDir, ...
    vLine, ex);

for i=numList
    [w1, ignore1, ignore2] = helper([ exDir1 num2str(i) exDir2 ], expDir, ...
        vLine, ex);
    s = reshape(w1, size(w1, 2)/ vLine, vLine);
    F = plotASample(s(markIndex, :), size(w1, 2)/ vLine, vLine, ...
        markList(markIndex));
    pdfFile = [expDir '/' exDir1 num2str(i) exDir2 '/results/sorted.pdf']
    saveAsPdf(F, pdfFile);
end
% l = dir(expDir);

% for i=3:size(l,1)
%     % Matrix file
%     [w, markIndex, markList] = helper(l(i), expDir, vLine, ex);
%     [w1, ignore1, ignore2] = helper(l(i), othExpDir, vLine, ex);
%     s = reshape(w1, size(w1, 2)/ vLine, vLine);
%     F = plotASample(s(markIndex, :), size(w1, 2)/ vLine, vLine, ...
%         markList(markIndex));
%     pdfFile = [othExpDir '/' l(i).name '/HebbPlot' l(i).name '.pdf'];
%     saveAsPdf(F, pdfFile);
%     close all;
% end

close all;

end


function [w, markIndex, markList] = helper(l, expDir, vLine, ex)
% Matrix file
matFile = [expDir '/' l '/results/matrix.txt'];

disp(matFile);

if(exist(matFile, 'file') == 2)
    disp(['Processing ' expDir '/' l '...']);
    % Read marks
    markFile = [expDir '/' l '/results/marks.txt'];
    id = fopen(markFile);
    markList = textscan(id, '%s', 'Delimiter' , ',');
    markList = reshape(markList{1}, 1, size(markList{1},1));
    markList = strrep(markList, '''', '');
    fclose(id);
    
    % Pdf file
    pdfFile = [expDir '/' l '/results/hebbPlotRearranged.pdf'];
    histFile = [expDir '/' l '/results/HistRearranged.pdf'];
    
    % This calculates the learning rate
    matrix = dlmread (matFile);
    [numRows,~] = size(matrix);
    learningRate = ((10000 / numRows) * 0.001);
    
    
    % Make plot
    [w, markIndex] = hebbPlot(matFile, vLine, learningRate, markList, ex, ...
        pdfFile, histFile);
    wFile = [expDir '/' l '/results/w.mat'];
    save(wFile, 'w');
    %close all;
else
    disp(matFile);
end
end