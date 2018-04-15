function markList = readMarks( markFile , removeSingleQuotes)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
% This function will read the marks from markFile and put it
% into a cell array. Will remove single quotes if specified.
% It will remove single quotes if they surround the marks in the file.

id = fopen(markFile);
markList = textscan(id, '%s', 'Delimiter' , ',');
markList = reshape(markList{1}, 1, size(markList{1},1));
if  removeSingleQuotes == 1
markList = strrep(markList, '''', '');
end

fclose(id);
end

