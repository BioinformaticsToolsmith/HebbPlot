

function s = dotSim( a, b )
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%SIM Calculates the similarity of two vectors using the inner product.

% Avoid dividing by zero
if(norm(a) == 0)
    aNorm = a;
else
    aNorm = a / norm(a);
end

% Avoid dividing by zero
if(norm(b) == 0)
    bNorm = b;
else
    bNorm = b / norm(b);
end

s = dot(bNorm , aNorm);
end

