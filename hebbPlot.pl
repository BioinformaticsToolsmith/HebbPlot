#! /usr/bin/perl -w
#
# Author: Alfredo Velasco
# Modified by Hani Z. Girgis
# The Bioinformatics Toolsmith, The University of Tulsa
# Logo designed by http://patorjk.com/software/taag/#p=display&h=0&v=3&f=Ogre&t=HebbPlot
#
#
# This will run everything
# WARNING THAT YOU SHOULD REALLY READ
# BEFORE EVEN THINKING ABOUT RUNNING THIS SCRIPT:
# THIS WILL OVERRIDE DATA UNDER THE DESTINATION DIRECTORY
#
my $installDir = '/Users/alfredo/Projects/Histone/Perl/HebbPlotPublic';

use File::Basename qw(dirname);
use Cwd qw(abs_path);
use Perl::epiUtils dirname( dirname abs_path $0) . '/Perl';
use Perl::sampleRegions dirname( dirname abs_path $0) . '/Perl';
use Perl::calculateOverlap dirname( dirname abs_path $0) . '/Perl';

#use Perl::epiUtils;
#use Perl::sampleRegions;
#use Perl::calculateOverlap;

use File::Path qw( rmtree );


my $logo = '		    __  __     __    __    ____  __      __ 
/ / / /__  / /_  / /_  / __ \/ /___  / /_
/ /_/ / _ \/ __ \/ __ \/ /_/ / / __ \/ __/
/ __  /  __/ /_/ / /_/ / ____/ / /_/ / /_  
/_/ /_/\___/_.___/_.___/_/   /_/\____/\__/ 
';

print $logo;


my $usage = "\n"
    . "Welcome to HebbPlot!\n"
    . "HebbPlot learns and displays a chromatine signature representing thousands of regions.\n"
    . "It is developed by Alfredo Velasco II and Hani Z. Girgis\n"
    . "The Bioinformatics Toolsmith Laboratory\n"
. "The University of Tulsa\n"
    . "Oklahoma, USA" . "\n" . "\n"
    . "Usage: $0 region_file "
    . "epigenome_directory "
    . "destination_directory "
    . "can_sort_mark_regions "
    . "number_of_regions_per_chromosome "
    . "expantion_ratio "
    . "number_of_vertical_lines "
    . "number_of_vertical_lines_to_exclude_at_peripheries "
    . "[colormap]\n" . "\n" . "\n"
    . "\tregion_file includes a list of regions of interest\n"
    . "\tepigenome_directory inlcudes chromatin mark files, each of which ends with .tagAlign\n"
    . "\t\tFor example: E031-H3K27me3.tagAlign and E031-H3K36me3.tagAlign\n"
    . "\tdestination_directory will include intermediate files and the HebbPlot in PDF form\n"
    . "\tcan_sort_mark_regions if the mark regions need to be sorted 1, otherwise 0\n"
    . "\tcan_merge_mark_regions if the mark regions need to be merged 1, otherwise 0\n"
    . "\tnumber_of_regions_per_chromosome is the number of samples taken per chromosome, e.g. 500\n"
    . "\texpansion_ratio is the ratio by which the regions will be expanded. It cannot be less than 0\n"
    . "\t\tFor example: a ratio of 0.1 will treat region:100-200 as region:90-210\n"
    . "\tnumber_of_vertical_lines will determine how many uniform samples will be taken from a region, e.g. 41\n"
    . "\tnumber_of_vertical_lines_to_exclude_at_peripheries is the number of samples to be excluded while sorting mark rows, e.g. 5\n"
    . "\tcolormap is the colormap to use for the final plot, e.g. gray\n"
    . "\tcolormap is an optional parameter and will be set to jet if ommited\n"
    . "\tThese are all valid colormaps: autumn bone colorcube cool copper flag gray grey hot hsv jet lines parula pink prism spring summer white winter\n"
    . "\n"
    . "Example\n"
    . "./hebbPlot.pl ./Enhancers/regions_enh_E066.bed ./Epigenomes/E066/ ./ExpEnhancers/ 1 1 500 0.1 41 5"
    . "\n" . "\n"
    . "Note: all files under the destination directory will be deleted before running the program\n"
    . "\n" . "\n";

die $usage unless ($#ARGV == 8 or $#ARGV == 9);

my (
    $regionFile, $epiDir,    $dstDir,   $canSort, $canMerge,
    $numSamples, $expansion, $numVLine, $numELine, $colorMap
    ) = @ARGV;

validate();

my $sortedEpiDir = "$dstDir/sortedEpigenome";
my $mergedEpiDir = "$dstDir/mergedEpigenome";
my $resultDir    = "$dstDir/results";
my $overlapDir   = "$dstDir/overlap";
my $sampleDir    = "$dstDir/sample";
my $sampleFile   = undef;
my $expandedFile = undef;

makeHebbPlot();

#
# Validate the user's paramters
#
sub validate {
    if ( !-e $regionFile ) {
	die "\n$regionFile does not exist. Please provide a valide region file\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    else {
	$regionFile = File::Spec->rel2abs($regionFile);
    }
    
    if ( !-d $epiDir ) {
	die "\n$epiDir does not exist. "
	    . "Please provide a valide epigenome directory\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    else {
	$epiDir = File::Spec->rel2abs($epiDir);
    }
    
    $dstDir = File::Spec->rel2abs($dstDir);
    
    if ( !( $canSort == 1 or $canSort == 0 ) ) {
	die "\ncan_sort_mark_regions must be 0 or 1\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if( !($canMerge == 1 or $canMerge == 0)){
	die "\ncan_merge_mark_regions must be 0 or 1\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if ( $numSamples < 1 ) {
	die "\nnumber_of_regions_per_chromosome must be at least one\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if ( $expansion < 0 ) {
	die "\nexpansion_ratio must be at least zero\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if ( $numVLine < 1 ) {
	die "\nnumber_of_vertical_lines must be at least one\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if ( $numELine < 0 ) {
	die "\nnumber_of_vertical_lines_to_exclude_at_peripheries "
	    . "cannot be negative\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    if ( $numVLine - 2 * $numELine <= 0 ) {
	die "\nnumber_of_vertical_lines_to_exclude_at_peripheries " /
	    "is too big for number_of_vertical_lines\n"
	    . "Please run the program without parameters for usage information\n"
	    . "\n";
    }
    
    my %colorMapList = ('parula' => 1, 'jet' => 1, 'hsv' => 1, 'hot' => 1,
			'cool' => 1, 'spring' => 1, 'summer' => 1,
			'autumn' => 1, 'winter' => 1, 'gray' => 1, 'grey' => 1, 'bone' => 1, 'copper' => 1, 'pink' => 1, 'lines' => 1,
			'colorcube' => 1, 'prism' => 1, 'flag' => 1, 'white' => 1);

    if (!defined($colorMap)){
	   $colorMap = 'jet';
    } else{
        $colorMap = lc($colorMap);
        if($colorMap eq 'grey'){
            $colorMap = 'gray'; # For our English friends
        }

        if(!exists($colorMapList{$colorMap})){
	    my @colorMapKeys = keys %colorMapList;
	    @colorMapKeys = sort(@colorMapKeys);
	    die("\n$colorMap is not a valid color map!\n"
		. "Valid color map options are as follows:\n"
		. "@colorMapKeys\n");
	}
    }
    
}

# This is the subroutine that runs all of the operations that make up hebbPlot.
# This includes preparing directories, taking samples, expanding, etc.
# It also prints the steps performed.
sub makeHebbPlot {
    print "HebbPlot will perform the following steps: \n";
    print "\t1) Prepare directories\n";
    print "\t2) Sort the input epigenome if requested by the user\n";
    print "\t3) Merge overlapping segments of the epigenome\n";
    print "\t4) Select a number of regions as specified by the user\n";
    print "\t5) Expand the input regions according to the ratio";
    print "\t6) Calculate the overlap between the marks and regions\n";
    print "\t7) Select uniform samples from each region\n";
    print "\t8) Generate the Hebb Plot, Quality Histogram, and Weight Vector\n";
    
    print "1) Preparing directories ...\n";
    prepareDirectories();
    
    print "2) Sorting if requested ...\n";
    if ( $canSort == 1 ) {
	sortEpigenome( $epiDir, $sortedEpiDir );
    }
    else {
	rmtree($sortedEpiDir);
	$sortedEpiDir = $epiDir;
    }
    
    print "3) Merging the epigenome ...\n";
    if($canMerge == 1){
	mergeDir( $sortedEpiDir, $mergedEpiDir );
    } else {
	rmtree($mergedEpiDir);
	$mergedEpiDir = $sortedEpiDir;
    }
    
    print "4) Selecting regions ...\n";
    makeSample();
    
    print "5) Expanding regions if requested ...\n";
    expandSampleFile();
    
    print "6) Calculating overlap between marks and regions ...\n";
    bedDir( $mergedEpiDir, $expandedFile, $overlapDir );
    
    print "7) Selecting uniform samples from the regions ...\n";
    initIntersection( $expandedFile, $overlapDir, $numVLine, $resultDir );
    
    print "8) Generating the Hebb Plot and Quality Histogram ...\n";
    runMatlab( $resultDir, $numVLine, $numELine );
}

# This subrouting deletes the destination directory
# and makes the subdirectories needed for outputs
sub prepareDirectories {
    # Delete the destination directory.
    if ( -d $dstDir ) {
	rmtree($dstDir) or die "Cannot remove $dstDir:$!\n";
    }
    mkdir($dstDir) or die "Cannot make $dstDir: $!\n";
    
    mkdir($sortedEpiDir) or die "Cannot make $sortedEpiDir: $!\n";
    mkdir($mergedEpiDir) or die "Cannot make $mergedEpiDir: $!\n";
    mkdir($resultDir)    or die "Cannot make $resultDir: $!\n";
    mkdir($overlapDir)   or die "Cannot make $overlapDir: $!\n";
    mkdir($sampleDir)    or die "Cannot make $sampleDir: $!\n";
}

# This subroutine will make the sample file
sub makeSample {
    my @regionSplit = split( '/', $regionFile );
    $sampleFile = "$sampleDir/$regionSplit[-1]";
    print "$regionFile $sampleFile\n";
    sampleFile( $regionFile, $sampleFile, $numSamples );
}

# This subroutine will expand the sample file
sub expandSampleFile {
    my @regionSplit = split( '/', $regionFile );
    $expandedFile = "$sampleDir/expanded_$regionSplit[-1]";
    if ( $expansion > 0.0 ) {
	expandRegionFile( $sampleFile, $expansion, $expandedFile );
    }
    else {
	$expandedFile = $sampleFile;
    }
}

# This subroutine calls Matlab to generate the HebbPlot
sub runMatlab {
    my $command = "matlab -nodisplay -nosplash -r "
	. "'cd $installDir/Matlab; generateHebbPlot('\\''$resultDir'\\'', $numVLine, $numELine, $colorMap); exit;'";
    print "$command\n";
    if ( system($command) ) {
	die "Could not execute $command: $!\n";
    }
}


