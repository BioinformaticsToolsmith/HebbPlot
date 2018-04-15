#! /usr/bin/perl -w

#
# Author: Alfredo Velasco
# The Bioinformatics Toolsmith, The University of Tulsa
#

# This utility module will download epigenetics files given a link and a file with
# the

#package Foo;
use Cwd;

require Exporter;
@ISA = qw(Exporter);
@EXPORT =
qw(mergeCoords min max mergeDir);

# This will eliminate overlaps in a given $inputFile and print out the results into an $outputFile
# The arguments are taken in as ($inputFile, $outputFile)
#
#      It turns this
#
#      ------------     -    ----            -----------      ------     ----------------- -------------
#      ----      ------     -------     - ----       ------   ------------     ------     ------------
#
#      into this
#
#      ---------------- -   -------     - -----------------   ------------------------------------------
sub mergeCoords {

	( my $repeatsInFile, my $repeatsOutFile ) = @_;

	open( IN, "<", $repeatsInFile )
	or die "Cannot open $repeatsInFile: $!\n";

	open( OUT, ">", "$repeatsOutFile" )
	or die "Cannot create $repeatsOutFile: $!\n";

	$_ = <IN>;

	my @line = split;
	while (<IN>) {
		my @temp = split;

		if ( @temp == 0 ) { next; }

		if ( $line[0] eq $temp[0]
			&& isOverlapping( $line[1], $line[2], $temp[1], $temp[2] ) )
		{
			$line[1] = min( $line[1], $temp[1] );
			$line[2] = max( $line[2], $temp[2] );
		}
		else {

			my $result = join( "\t", @line );
			print OUT "$result\n";
			@line = @temp;

		}

	}

	my $result = join( "\t", @line );
	print OUT "$result\n";

	$result = join( "\t", @temp );
	print OUT "$result\n";

	close(IN);
	close(OUT);

}


# returns whether or not the epigenomes are overlapping
# based on their starts and ends
# arguments given as ($start1, $end1, $start2, $end2)
sub isOverlapping {
	my ( $s1, $e1, $s2, $e2 ) = @_;

	my $isStartWithin = ( ( $s2 >= $s1 ) && ( $s2 <= $e1 ) );
	my $isEndWithin   = ( ( $e2 >= $s1 ) && ( $e2 <= $e1 ) );
	my $isIncluding   = ( ( $s2 >= $s1 ) && ( $e2 <= $e1 ) );
	my $isIncluded    = ( ( $s1 >= $s2 ) && ( $e1 <= $e2 ) );
	my $isAdjacent = ( $e1 + 1 == $s2 ) || ( $e2 + 1 == $s1 );

	return ( $isStartWithin
		|| $isEndWithin
		|| $isIncluding
		|| $isIncluded
		|| $isAdjacent );
}

# Given a directory, it will perform mergeCoords() on every file
# Arguments given as ($parentDir, $destDir)
# 
sub mergeDir {
	my ( $fol, $dstDir ) = @_;
	my @files = glob( $fol . '/*.tagAlign' );
	foreach my $file (@files) {
		my @splitFile = split( '/', $file );
		my $nickName = $splitFile[-1];

		if ( $nickName =~ m/.*\d.*tagAlign$/ ) {

			if ( $nickName =~ m/(.*)(\.tagAlign)$/ ) {
				$outFile = "$dstDir/$1M$2";
			}

			print "Merging $file $outFile\n";

			mergeCoords( $file, $outFile );
		}
	}

}

# This runs the unix sort command on epigenomes found
# in a given directory and output the results into another directory.
# Arguments give as (epigenomes directory, output directory)
sub sortEpigenome {
	my ( $epiDir, $sortedEpiDir ) = @_;
	my @markList = glob( $epiDir . "/*.tagAlign" );
	foreach my $mark (@markList) {
		print "\tSorting $mark ...\n";
		my @splitMark = split( '/', $mark );
		my $sortCom =
		"sort -k1,1 -k2,2n $mark > $sortedEpiDir/$splitMark[-1]\n";
		if ( system($sortCom) ) {
			die "Could not execute $sortCom: $!\n";
		}
	}
}

# Simple minimum function
# Arguments give as ($num1, $num2)
sub min {
	my ( $a, $b ) = @_;
	if   ( $a > $b ) { return $b; }
	else             { return $a; }
}

# Simple maximum function
# Arguments give as ($num1, $num2)
sub max {
	my ( $a, $b ) = @_;
	if   ( $a < $b ) { return $b; }
	else             { return $a; }
}

#1;
