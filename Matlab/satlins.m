
function a = satlins( b )
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%

%SATLINS calculate the symmetric saturating linear function

a = b;
a(b < -1) = -1;
a(b > 1) = 1;
end

