function [S, markIndex] = sortASampleHierarchical(m, r, c, e)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%

%SORTASAMPLEHIERARCHICAL sort a matrix using hierarchical clustering
% m is the matrix
% r is the number of rows
% c is the number of columns
% e is the number of columns to EXCLUDE on each side
%

H = reshape(m, r, c);
if r > 1
    A = H;
    if e > 0
        A = H(:, e:c-e);
    end
    
    cityBlock = 'cityblock';
    Z = linkage(A, 'weighted', cityBlock);
    D = pdist(A, cityBlock);
    markIndex = optimalleaforder(Z, D);
    O = ones(1,size(A,2));
    simToFirst = dotSim(O, A(markIndex(1), :));
    simToLast = dotSim(O, A(markIndex(end), :));
    
    if simToFirst >= simToLast
        markIndex = markIndex(end:-1:1);
    end
    S = H(markIndex, :);
else
    S = H;
    markIndex = [1];
end

end

