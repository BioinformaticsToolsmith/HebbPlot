#! /usr/bin/perl -w

#
# Author: Alfredo Velasco
# The Bioinformatics Toolsmith, The University of Tulsa
#
# This program samples region files

require Exporter;
@ISA    = qw(Exporter);
@EXPORT = qw(sampleFile);

use strict;
use List::Util qw(shuffle);

my $inDir;
my $outDir;
my $numSamples;


# Given a region file, output file, and number of samples
# this subroutine will read the region file and gather the regions
# by chromosome. It will select samples given by the number of samples
# parameter and output the samples into the output file
sub sampleFile {
	my ( $inFile, $outFile, $numSamples ) = @_;
	if ($numSamples <= 0 or ! ($numSamples =~ /^\d+\z/)){
		die ("$numSamples is not a valid number of samples! The number of samples must be a positive integer\n");
	}

	#print "I am this one\n";

	open( IN,  "<$inFile" )  or die "Cannot open $inFile: $!\n";
	open( OUT, ">$outFile" ) or die "Cannot open $outFile: $!\n";
	my @chromosomeList = ();
	while (<IN>) {
		chomp;
		my $readLine = $_;

		my @initialLine = split("\t",$readLine);

#Modification made by Zachary Reyes on 07/18/17.
# The following line was added to make sure that the sampleFile subroutine was able to handle bed files that had
# extra information on each line.
# It does this simply by grabbing the first 3 space-delimited elements in each bed file.
# The resulting bed file produced by this tool will only include the first 3 elements in each line.
		my @line = @initialLine[ 0, 1, 2 ];
#		print join("\t", @line) . "\n";
		if ( $#chromosomeList == -1 and not ($readLine eq "") ) {
			push @chromosomeList, \@line;
		}
		elsif ( @line == 0 or $readLine eq "") {
			# print "Empty line!\n";

			# When the line is empty
		}
		elsif ( $line[0] eq $chromosomeList[-1]->[0] ) {
			push @chromosomeList, \@line;
		}
		else {
			my $interval = int( $#chromosomeList / $numSamples );
			if ( $interval == 0 ) {
				$interval = 1;
			}
			for ( my $i = 0 ; $i <= $#chromosomeList ; $i = $i + $interval ) {
				my @result = @{ $chromosomeList[$i] };
				print OUT join "\t", @result;
				print OUT "\n";
			}
			@chromosomeList = ();
			push @chromosomeList, \@line;
		}
	}
	my $interval = int( $#chromosomeList / $numSamples );
	if ( $interval == 0 ) {
		$interval = 1;
	}
	for ( my $i = 0 ; $i <= $#chromosomeList ; $i = $i + $interval ) {
		my @result = @{ $chromosomeList[$i] };
		print OUT join "\t", @result;
		print OUT "\n";
	}
	@chromosomeList = ();
	close(OUT);
	close(IN);
}

1;
