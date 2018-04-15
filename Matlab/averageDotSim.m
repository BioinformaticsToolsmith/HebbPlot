
function averageDotSim(dir1, dir2, minimum)
% AVERAGEDOTSIM
% Authors: Hani Z Girgis and Alfredo Velasco
% The Bioinformatic Toolsmith Laboratory
% The University of Tulsa
% This function takes in two directories
% that represent a tissue.
% It will go through the w.mat files (the weights
% associated with the hebbPlots) and run 
% dotSim on corresponding marks in the same tissue
% under the two directories.
% The minimum will remove marks with counts less than it
% This will output the results to the console
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
l = dir(dir1);

tableAll = containers.Map({'All Marks'}, {0});
tableDotSim = containers.Map({'Dot Sim'}, {0.0});

for i=3:size(l,1)
    % Read marks
    name = l(i).name;
    disp([dir1 name '/results/marks.txt']);
    markList = readMarks([dir1 name '/results/marks.txt'], 0);
    
    % Count the number of occurrances of all marks
    for j=1:size(markList, 2)
        key = char(markList(j));
        
        if ~isKey(tableAll, key)
            tableAll(key) = 1;
        else
            tableAll(key) = tableAll(key) + 1;
        end
    end
    
    % Update the dotSim table
    mHigh = importdata([dir1 name '/results/w.mat']);
    disp([dir2 name '/results/w.mat']);
    disp([dir2 name '/results/w.mat']);
    mLow = importdata([dir2 name '/results/w.mat']);
    c = 41;
    r = size(mHigh,2)/c;
    high = reshape(mHigh, r, c);
    low = reshape(mLow, r, c);
    for j = 1:r
        key = char(markList(j));
        if ~isKey(tableDotSim, key)
            if strcmp(key, 'H4K12ac')
                disp( name );
            end
            tableDotSim(key) = dotSim(high(j,:), low(j,:));
        else
            tableDotSim(key) = tableDotSim(key) + dotSim(high(j,:), low(j,:));
        end
    end
    
end

% After this loop, the dot sim table contains the average not the total
keyList = keys(tableDotSim);
for j=1:size(keyList,2)
    key = keyList{j};
    if(~strcmp(key, 'Dot Sim'))
        tableDotSim(key) = tableDotSim(key) / tableAll(key);
    end
end

% Sort the dot sim table
valueList = cell2mat(values(tableDotSim));
[s, i] = sort(valueList);
sortedKeyList = keyList(i);
disp(['Mark' '    ' 'Count' '    ' 'DotSim']);
for j=1:size(sortedKeyList,2)
    key = sortedKeyList{j};
    if(~strcmp(key, 'Dot Sim'))
        if tableAll(key) >= minimum
            disp([key '    ' num2str(tableAll(key)) '    ' num2str(tableDotSim(key))]);
        end
    end
end

% allMarks = keys(tableAll);
% sortedMarks = zeros(size(allMarks, 2), 1);
% 
% for j=1:size(allMarks, 2)
%     disp( [char(allMarks(j)) '     ' tableAll(char(allMarks(j)))]);
% end
% % pause;
% for j=1:size(allMarks, 2)
%     key = char(allMarks(j));
%     disp(key);
%     sortedMarks(j:1) = key;
%     sortedMarks(j:2) = num2str(tableDotSim(key) / tableAll(key));
% end
% disp(sortedMarks);
% 
% for j=1:size(allMarks, 2)
%     key = char(allMarks(j));
%     if isKey(tableDotSim, key)
%         disp([key  '    ' ...
%             num2str(tableAll(key))  '    ' ...
%             num2str(tableDotSim(key))  '    ' ...
%             num2str(tableDotSim(key) / tableAll(key) )]);
%     end
% end

end