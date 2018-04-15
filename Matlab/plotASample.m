

function h = plotASample( m, r, c, markLabels, colormapIn, varargin)
%
% Author: Hani Zakaria Girgis, PhD
% The Bioinformatics Toolsmith Laboratory
% The University of Tulsa
%heb
%PLOTASAMPLE convert a vector to a matrix. Then, it makes an image of it.

s = reshape(m, r, c);

% Plot the unsorted matrix
if nargin > 6
    h = figure; imagesc(s, varargin{1}); colormap(colormapIn);
else
    h = figure; imagesc(s, [-1 1]); colormap(colormapIn);
end

colorbar;

%Set the tick labels
part1 = floor(c/2):-1:1;
part2 = 1:floor(c/2);
xTick = 1:c;

xTickLabels = [part1, 0, part2];
q = floor(c/4);
center = round(c/2);
set(gca, 'XTick', xTick([1 center-q center center+q c]), 'XTickLabel', ...
    xTickLabels([1 center-q center center+q c]), ...
    'TickLabelInterpreter', 'latex');

set(gca, 'YTick', 1:r, 'YTickLabel', markLabels);
set(gca, 'FontName', 'Times');


xlabel('Uniform Samples', 'interpreter', 'latex', 'FontName', 'Times');
ylabel('Marks', 'interpreter', 'latex', 'FontName', 'Times');

% Set the font size of all text in the figure
set(gca, 'FontSize', 14);
end