
function comparePosNegPromters(nDir, pDir, minimum)
% COMPAREPOSNEGPROMTERS
% Authors: Hani Girgis and Alfredo Velasco
% The Bioinformatic Toolsmith Laboratory
% The University of Tulsa
% This function takes in two directories
% that represent a tissue and a minimum.
% It will go through the w.mat files (the weights
% associated with the hebbPlots) and calculates
% the marks that have directional preferences
% The minimum will remove marks with counts less than it
% Directorys must have the following structure
% Dir
%  -> Epigenome001
%    -> results
%      -> files
%  -> Epigenome002
%    -> results
%      -> files
%  ...
% Directories must also have the same names for subdirectories
% Warning: The directories must have slashes at the end
l = dir(pDir);
tableOpposite = containers.Map({'Directional Marks'},{0});
tableAll = containers.Map({'All Marks'}, {0});

for i=3:size(l,1)
    name = l(i).name;
    disp([pDir name '/results/marks.txt']);
    markList = readMarks([pDir name '/results/marks.txt'], 0);
    
    % Count the number of occurrances of all marks
    for j=1:size(markList, 2)
        key = char(markList(j));
        if ~isKey(tableAll, key)
            tableAll(key) = 1;
        else
            tableAll(key) = tableAll(key) + 1;
        end
    end
    
    mP = importdata([pDir name '/results/w.mat']);
    mN = importdata([pDir name '/results/w.mat']);
    
    
    c = 91;
    cThird = 31;
    r  = size(mP,2)/c;
    h  = reshape(mP, r, c);
    hL = reshape(h(:,1:31), 1, cThird  * r);
    hR = reshape(h(:,61:91), 1, cThird  * r);
    
    [s, i] = sortDifferentMarks( hL, hR, cThird);
    
    oppositeMarks = markList(i(s < 0));
    oppositeMarksNum = size(oppositeMarks, 2);
    
    % Count the number of marks showing directional preferences
    for j =1:oppositeMarksNum
        key = char(oppositeMarks(j));
        if ~isKey(tableOpposite, key)
            tableOpposite(key) = 1;
        else
            tableOpposite(key) = tableOpposite(key) + 1;
        end
    end
end

% Print: Mark Occurrances Directional
disp(['Mark ' 'Known ' 'Directional ' 'Percentage']);
allMarks = keys(tableAll);
for j=1:size(allMarks, 2)
    key = char(allMarks(j));
    if isKey(tableOpposite, key)
        if tableAll(key) >= minimum
            disp([key  ',' ...
                num2str(tableAll(key)) ',' ...
                num2str(tableOpposite(key)) ',' ...
                num2str(round(100 * tableOpposite(key)/tableAll(key)))]);
        end
    end
end

end
