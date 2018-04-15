

function w = trainOutstar( p0, p, alpha )
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%
%trainOutstar This function learns a prototype
%   It uses the outstar hebb rule

w0 = ones(1, size(p0,2));
w = zeros(1, size(p0,2));
% convergenceList = zeros(1, size(p0, 1));

for i = 1:size(p0, 1)
% wOld = w;
a = satlins((p0(i, :) .* w0) + (p(i,:) .* w) );
w = w + alpha * (a - w) .* p(i, :);
% convergenceList(i) = dotSim(w, wOld);
end
% plot(convergenceList(1:500));
% print ('/Users/alfredo/Writing/HebbPlot/Supplementary/Convergence/Promoters_E003.pdf', '-dpdf');
% error('nothing')
end

