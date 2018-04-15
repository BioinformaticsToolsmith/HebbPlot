#! /usr/bin/perl -w

#
# Author: Alfredo Velasco
# The Bioinformatics Toolsmith, The University of Tulsa
#

# This utility module holds the processes needed to compute the regional overlaps
# It prints 3 files: one for the matrix, one for the regions(the matrix keys), and one for the marks(E001-H3H327B.tagAlign, E001-H3H3lkjsdf.tagAlign, etc.)
# They all must be in the same order

use Cwd;

use Perl::epiUtils dirname( dirname abs_path $0) . '/Perl';

use List::MoreUtils qw(uniq);
use List::Flatten;
use Math::BigFloat ':constant';
use Data::Dumper;
use File::Path qw( rmtree );

use strict;

require Exporter;
my @ISA    = qw(Exporter);
my @EXPORT = qw(downloadRegions bedDir);

my %matrix         = ();
my $regFile        = undef;    # Is only needed to initialize the matrix
my $epiDir         = undef;
my $numVLines      = undef;
my $dstDir         = undef;
my $marksProcessed = 0;

# This must be called first. It will initialize
# ( region file, epigenome directory, number of vertical lines, output directory)
sub initIntersection {

	my ( $regFileLoc, $epiDirLoc, $numVLinesLoc, $outputDir ) = @_;
	if ( !-e $regFileLoc ) {
		die("$regFileLoc does not exist!");
	}
	if ( !-d $epiDirLoc ) {
		die("$epiDirLoc does not exist!");
	}
	$regFile        = $regFileLoc;
	$epiDir         = $epiDirLoc;
	$numVLines      = $numVLinesLoc;
	$dstDir         = $outputDir;
	$marksProcessed = 0;
	initMatrix();
	processEpigenome();
	printMatrix();

}

# This reads the region file
# and stuffs the regions as keys
# into a hash table 
sub initMatrix {
	if ( defined $regFile ) {
		;
		open( IN, $regFile ) or die("Could not open $regFile: $!\n");
		%matrix = ();
		while (<IN>) {
			my @t = split;
			if ( @t != 0 ) {
				$matrix{ $t[0] . "-" . $t[1] . "-" . $t[2] } = [];

			}
		}
		close(IN);
	}
	else {
		die( "The region file is undefined.\nPlease call init(regFile, epiDir)."
			);
	}
}

# Given a folder like E001 (which should have been passed when init was called),
# it gets all of the E001H3H3-23S.tagAlign files and calls processMarkFile on all of them.
# Arguments given as (epigenomeDirectory)
sub processEpigenome {
	my @files = glob( $epiDir . '/*.tagAlign' );
	if ( $#files == -1 ) {
		die( "Nothing was globbed from $epiDir" . '/*.tagAlign' );
	}

	my $markFile = $dstDir . "/marks.txt";
	open( MARKS, "> $markFile" ) or die "Cannot open $markFile: $!";
	my @output = ();

	foreach my $file (@files) {

		print "Generating the matrix from $file ...\n";
		if ( $file =~ m/.*\/E\d\d\d-(.*)\.tagAlign/ ) {
			# print "\n\n$file $1 is file\n\n";
			push( @output, $1 );
		}

		processMarkFile($file);
	}

	print MARKS join ',', @output;

	close(MARKS);

}

# Given a mark file, it puts the information in %regions
sub processMarkFile {
	my ($markFile) = @_;
	$marksProcessed++;

	my ( $region, @regionList );
	open( IN, '<', $markFile ) or die "Cannot open $markFile: $!\n";

	my $intersect      = ();
	my %intersectTable = ();

	while (<IN>) {
		my @temp   = split;
		# Works with original bedtools files
		my $region = $temp[6] . "-" . $temp[7] . "-" . $temp[8];
		my $mark   = [ $temp[0], $temp[1], $temp[2] ];
		if ( !exists $intersectTable{$region} ) {
			$intersectTable{$region} = [];
		}

		push $intersectTable{$region}, $mark;

	}
	close(IN);

	# We will push the line intersection row into the hash
	while ( my ( $region, $overLapList ) = each %intersectTable ) {
		my $result = lineIntersection( $region, $overLapList );

		if ( exists $matrix{$region} ) {
			push $matrix{$region}, $result;
		}
		else {
			die "Region $region is not in the matrix\n";
		}
	}

	# Any region that didn't intersect with anything
	# will be represented as a row of -1s
	while ( my ( $key, $value ) = each(%matrix) ) {
		if ( $#$value != $marksProcessed - 1 ) {
			my @result = fillEmpty();
			push @{ $matrix{$key} }, \@result;
		}
	}

	# Post-condition: make sure that the number of columns is the same as the marks
	while ( my ( $key, $value ) = each(%matrix) ) {
		if ( $#$value != $marksProcessed - 1 ) {
			print Dumper($value) . "\n$key\n";
			die("Too many items!");
		}
	}
}

# This will print the matrix column-wise
#
#	Example:
#	1	2	3
#	4	5	6
#	7	8	9
#
#   will be printed as
#
#	1	4	7	2	5	8	3	6	9
#
#
sub printMatrix {

	my ( $file1, $file2 ) = ( "$dstDir/matrix.txt", "$dstDir/regions.bed" );
	open( OUTM, ">$file1" );
	open( OUTR, ">$file2" );

	while ( my ( $key, $value ) = each %matrix ) {
		print OUTR ("$key\n");
		my @result = flatten2DArray($value);
		if ( $#result + 1 != $numVLines * $marksProcessed ) {
			die("Number of items in matrix row is wrong! @result");
		}
		print OUTM "@result\n";
	}
	close OUTR;
	close OUTM;
}

# Given a reference to a 2D array, it will return them
# column-wise.
#
#	Example:
#	1	2	3
#	4	5	6
#	7	8	9
#
#   will be printed as
#
#	1	4	7	2	5	8	3	6	9
#
# Arguments given as ($arrayReference)
# This will return a new array or a reference to the array of results
sub flatten2DArray {
	my ($array) = @_;

	my @result = ();
	my $index  = 0;
	my $used   = 1;    # If no elements are printed out in this iteration,
	                   # that means we' ve gone through every item

	# in the array
	while ($used) {
		$used = 0;
		foreach my $rows (@$array) {
			if ( $index <= $#$rows ) {
				push @result, $$rows[$index];
				$used = 1;    # We've printed an item
			}
		}
		$index++;
	}
	return wantarray ? @result : \@result;
}

# Returns a list of -1's or a reference.
# This only works if the list given is empty
# Arguments given as ($arrayRef)
sub fillEmpty {
	my @result = ();

	for ( my $i = 0 ; $i < $numVLines ; $i++ ) {
		push @result, -1;
	}
	if ( $#result != $numVLines - 1 ) {
		die("-1's has is the wrong length");
	}

	#	print "$#$array $#result\n";
	return wantarray ? @result : \@result;
}


# Downloads the regions given the main link that holds the files, the input file with the names of the regions you want to download,
# and the directory where you want them stored
# Arguments taken in as ($parentLink, $inputFile, $dstDirectory)
sub downloadRegions {
	( my $parentLink, my $inputFile, my $parentDir ) = @_;
	if ( -d $parentDir ) { }
	else {
		mkdir($parentDir);
	}
	open( IN, $inputFile ) or die "Could not open $inputFile: $!\n";
	while (<IN>) {
		my $link = (split)[0];
		my $cmd  = "curl \"$parentLink$link\" -o $parentDir/$link";
		print "$cmd\n";
		if ( system($cmd) ) {
			die("Could not execute $cmd: $!");
		}
	}
	close(IN);
	decompress($parentDir);
}


# Given a region directory and an expansion rate
# this subroutine will read the region files
# in the region directory, then read the regions
# and output their expanded versions by
# performing the operation on them:
# length = end - start;
# start = start + length * expansion
# end = end - length * expansion
sub expandRegion {
	my ( $regionDir, $expansion ) = @_;
	my @regionList = glob( $regionDir . "*" );
	my @regDir = split( '/', $regionList[0] );
	pop @regDir;
	my $expRegDir = join( '/', @regDir ) . "X";
	@regDir = split( '/', $regionList[0] );
	if ( -d $expRegDir ) {
		rmtree($expRegDir);
	}
	mkdir($expRegDir);

	foreach my $reg (@regionList) {
		open( IN, '<' . $reg ) or die "Could not open $reg: $!\n";
		my @file = split( '/', $reg );
		open( OUT, '>' . "$expRegDir/$file[-1]" )
		or die "Could not open $reg: $!\n";
		print "Expanding $expRegDir/$file[-1]\n";
		while (<IN>) {
			my @temp     = split();
			my $interval = $temp[2] - $temp[1];
			my $delta    = $interval * $expansion;
			$temp[2] = int( $temp[2] + $delta );
			$temp[1] = int( $temp[1] - $delta );
			if ( $temp[1] < 0 ) {
				$temp[1] = 0;
			}
			print OUT join( "\t", @temp );
			print OUT "\n";
		}

		close(IN);
		close(OUT);

		#		last;
	}
	return $expRegDir;
}

# Given a region file, an expansion rate,
# and an output file, 
# this subroutine will read the regions
# and output their expanded versions by
# performing the operation on them:
# length = end - start;
# start = start + length * expansion
# end = end - length * expansion
# These expanded regions will be outputed
# to the output file
sub expandRegionFile {
	my ( $regionFile, $expansion, $outputFile ) = @_;

	open( IN, '<' . $regionFile ) or die "Could not open $regionFile: $!\n";
	open( OUT, '>' . "$outputFile" )
	or die "Could not open $outputFile: $!\n";
	while (<IN>) {
		# print "|$_|";
		my @temp     = split();
		if($#temp == -1){
			next;
		}
		my $interval = $temp[2] - $temp[1];
		my $delta    = $interval * $expansion;
		$temp[2] = int( $temp[2] + $delta );
		$temp[1] = int( $temp[1] - $delta );
		if ( $temp[1] < 0 ) {
			$temp[1] = 0;
		}
		print OUT join( "\t", @temp[ 0, 1, 2 ] );
		print OUT "\n";
	}

	close(IN);
	close(OUT);

}

# Given a directory and region.bed, it will execute the bedtools intersect command
# on the files in the directory and write the results into a file
# Arguments given as ($directory, $region.bed)
# Warning: It only works on files that end with M.tagAlign
sub bedDir {
	my ( $fol, $bed, $dstDir ) = @_;
	my @files = glob( $fol . "/*.tagAlign" );

	foreach my $file (@files) {
		if ( $file =~ m/.*\/(.*)M.tagAlign/ ) {
			my $outputFile = "$dstDir/$1.tagAlign";

			my $command =
			"\tbedtools intersect -a $file -b $bed -wb > $outputFile\n";
			print "$command";

			if ( system($command) ) {
				die "Cannot execute $command: $!";
			}
		}
	}

}

# This subroutine takes in a region and list refernce
# that holds areas of interest.
# Look the this example to see what I mean. Remember that {} and [] are references.
# | | | | | | | | | | | | | | | | | | | | | | | | | |
# v v v v v v v v v v v v v v v v v v v v v v v v v v
#
# my $testp1 = 'chr1:100-300' ; my $testp2 = [["chr1",100,120],["chr1",130,140],["chr1",160,210]];
#
# It takes in the region and uses the number of vertical lines to determine intersections.
# Example: if the region is chr:100-300 and the number of lines is 5, intersections
# will be tested out at points 100, 150, 200, 250, and 300
# Find the intersections
# Arguments give as ($regionOfInterestFile, $markListReference )
sub lineIntersection {
	my ( $region, $markListRef ) = @_;

	if ( !defined($numVLines) || $numVLines <= 0 ) {
		die(
			"The number of lines desired is invalid; call init() and pass acceptable arguments!\n"
			);
	}

	my ( $char, $start, $end ) = split( /[\-:]/, $region );

	my $dx = ( $end - $start + 1 ) / ( $numVLines - 1 );

	# Making a list of the indexes of the vertical lines
	my @vLineList = ($start);
	for ( my $h = 1 ; $h < $numVLines ; $h++ ) {
		if ( $h != $numVLines - 1 ) {
			push( @vLineList, $vLineList[-1] + $dx );
		}
		else {
			push( @vLineList, $end );
		}
	}
	for ( my $i = 0 ; $i <= $#vLineList ; $i++ ) {
		if ( $vLineList[$i] > $end ) {
			$vLineList[$i] = $end;
		}
	}

	my $index     = 0;
	my $lastIndex = $#$markListRef;

	my $vLineIndex = 0;
	my @result     = ();
	while ( $vLineIndex <= $#vLineList && $index <= $lastIndex ) {

   # This is what you need to know. It's the chr1, start, and end of the region.
   my $markRef = $markListRef->[$index];
   my ( $chr, $s, $e ) = @$markRef;

   if ( $vLineList[$vLineIndex] <= $e && $vLineList[$vLineIndex] >= $s ) {
   	push( @result, 1 );
   	$vLineIndex++;
   }
   elsif ( $vLineList[$vLineIndex] < $s ) {
   	push( @result, -1 );
   	$vLineIndex++;
   }
   elsif ( $vLineList[$vLineIndex] > $e ) {
   	$index++;
   }
}

while ( $index > $lastIndex && $vLineIndex <= $#vLineList ) {
	push( @result, -1 );
	$vLineIndex++;
}


if ( $#result == -1 ) {
	die("Result is empty1");
}

if ( $#result + 1 != $numVLines ) {
	die(    "Result length is "
		. ( $#result + 1 )
		. " instead of $numVLines!\n" );
}

return wantarray ? @result : \@result;
}

# This directory takes in two directories. One with the epigenomes, and one with region files in the .bed format
# This method will go through all of the region files and epigenomes and perform the intersection task on all of them.
# If you want to intersect a specific epigenome(E001, for example) use the bedDir method
# Arguments given as ($epigenomeDir, $bedtoolsDir, $destinationDir)
sub masterBedDir {
	my ( $fol, $bed, $dstDir, $expansion ) = @_;

	if ( !-d $dstDir ) {
		mkdir($dstDir) or die "Cannot make $dstDir: $!\n";
	}
	if ($expansion) {
		$bed = expandRegion( $bed, $expansion );
	}

	my @regionList = glob( $bed . "/*.bed" );

	for ( my $i = 0 ; $i <= $#regionList ; $i++ ) { 
		if ( $regionList[$i] =~ m/.*(E\d{3})_(\d+).*/ ) {
			my $epigenomeIn  = "$fol/$1";
			my $epigenomeOut = "$dstDir/$1_$2";
			if ( -d $epigenomeIn ) {

				bedDir( $epigenomeIn, $regionList[$i], $epigenomeOut );
			}
		}

	}
	return $bed;
}

1;
