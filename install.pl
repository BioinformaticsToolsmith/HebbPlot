#! /usr/bin/perl -w 
# Author: Alfredo Velasco
# This installs the hebbPlot.
# It does this by changing the current directory.
use File::Copy;
use CWD;
my $currentDir = `pwd`;
chomp($currentDir);

my $temp = 'temp.pl';
my $hebb = 'hebbPlot.pl';
my $isWrittenYet = 0;
open(IN, '<', $hebb) or die "Cannot open $hebb: $!\n";
open(OUT, '>', $temp) or die "Cannot open $temp: $!\n";

while(<IN>){
	if(! m/#/ and ! $isWrittenYet){
		print OUT 'my $installDir = ' . "'$currentDir'" . ';' . "\n";
		$isWrittenYet = 1;
	}
	if ( s/my .installDir = .*//){
		next;
	}
	
	print OUT;
}
close(OUT);
close(IN);

move($temp, $hebb);
chmod(0700, $hebb);

#unlink($temp);
