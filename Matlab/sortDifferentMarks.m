

function [s, i] = sortDifferentMarks( v1, v2, c, varargin)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%SortMarks this function sorts the marks starting from the most dissimilar
%to the most similar. It must be applied to signatures that have the same
%marks. For exmaple, compare the signatures of the repetitive and the 
%non-repetitive enhancers specific to the same tissue
%

r1 = size(v1,2) / c;
r2 = size(v2,2) / c;

m1 = reshape(v1, r1, c);
m2 = reshape(v2, r2, c);

v1RowIndex = 1:r1; 
v2RowIndex = 1:r2;

if nargin == 5
  v1RowIndex = varargin{1};
  v2RowIndex = varargin{2}; 
end

m1 = m1(v1RowIndex, :);
m2 = m2(v2RowIndex, :);

s = zeros(1,length(v1RowIndex));

for j=1:length(v1RowIndex)
    s(1,j) = dotSim(m1(j,:), m2(j,:));
end
[s, i] = sort(s);

end