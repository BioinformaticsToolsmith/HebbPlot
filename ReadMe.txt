# Authors:
# Alfredo Velasco: alv657@utulsa.edu
# Hani Girgis: hani-girgis@utulsa.edu


Requirements:

Perl (https://www.perl.org)
Bedtools (http://bedtools.readthedocs.io/en/latest/)
Matlab (https://www.mathworks.com)
Statistics and Machine Learning Toolbox (https://www.mathworks.com/products/statistics.html)
PDF Viewer (https://get.adobe.com/reader/)
UNIX/Linux/Mac 

To install, run the following command:
./install.pl

To learn about the paramaters of HebbPlot, run the following command:
./hebbPlot.pl

Note about sort and merge commands: if your epigenomes aren't sorted/merged, you can tell hebbPlot
to sort and merge them for you. After doing it once, you can reuse those same sorted/merged epigenomes.
Simply move the overlap/ directory to and specify this directory the next time you run hebbPlot again.


If you want to produce results similar to Table 1, as seen in the paper,
use comparePosNegProm.m in the Matlab directory.

If you want to produce results similar to Table 2 and 3, as seen in the paper,
use averageDotSim.m in the Matlab directory.

Usage: /Users/alfredo/Projects/Histone/Perl/HebbPlotPublic/hebbPlot.pl region_file epigenome_directory destination_directory can_sort_mark_regions number_of_regions_per_chromosome expantion_ratio number_of_vertical_lines number_of_vertical_lines_to_exclude_at_peripheries [colormap]


	region_file includes a list of regions of interest
	epigenome_directory inlcudes chromatin mark files, each of which ends with .tagAlign
		For example: E031-H3K27me3.tagAlign and E031-H3K36me3.tagAlign
	destination_directory will include intermediate files and the HebbPlot in PDF form
	can_sort_mark_regions if the mark regions need to be sorted 1, otherwise 0
	can_merge_mark_regions if the mark regions need to be merged 1, otherwise 0
	number_of_regions_per_chromosome is the number of samples taken per chromosome, e.g. 500
	expansion_ratio is the ratio by which the regions will be expanded. It cannot be less than 0
		For example: a ratio of 0.1 will treat region:100-200 as region:90-210
	number_of_vertical_lines will determine how many uniform samples will be taken from a region, e.g. 41
	number_of_vertical_lines_to_exclude_at_peripheries is the number of samples to be excluded while sorting mark rows, e.g. 5
	colormap is the colormap to use for the final plot, e.g. gray
	colormap is an optional parameter and will be set to jet if ommited
	These are all valid colormaps: autumn bone colorcube cool copper flag gray grey hot hsv jet lines parula pink prism spring summer white winter

Example
./hebbPlot.pl ./Enhancers/regions_enh_E066.bed ./Epigenomes/E066/ ./ExpEnhancers/ 1 1 500 0.1 41 5 jet

Note: all files under the destination directory will be deleted before running the program
