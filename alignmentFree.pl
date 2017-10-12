#!/usr/bin/env perl
use warnings;
use strict;
my $installDir = '/home/benjamin/Documents/Today/Dist';

#
# Author: Brian Luczak
# Author: Ben James
# Author: Hani Girgis
#

use strict;
use Cwd qw(abs_path);
use File::Path qw(rmtree);

die "Usage: $0 fasta_directory output_directory alignment_type" unless $#ARGV == 2;

my ( $faDir, $outDir, $altype ) = @ARGV;

my $delim       = '!';
my $matlabDir   = "$installDir/Matlab/code";
my $cppDir      = "$installDir/Cpp";
# use -l option to use local alignment for alignment-tool
my $aligner     = "$cppDir/Align/bin/alignment-tool";
my $histogramer = "$cppDir/Hist/bin/aligner";
my $histDir     = "$outDir/Hist";
my $alignDir    = "$outDir/Align";

drive();

sub drive {

	# Make full qualified directory names
	$faDir  = abs_path($faDir);
	$outDir = abs_path($outDir);
	getFreshDir($outDir);
	getFreshDir($histDir);
	getFreshDir($alignDir);

	runAlignHist();
	runMatlab();
}

sub runAlignHist {
	my $coverage = 0;
	if ($altype eq "global") {
	} elsif ($altype eq "local") {
		$aligner = "$aligner -l";
	} elsif ($altype eq "coverage") {
		$aligner = "$aligner -l";
		$coverage = 1;
	} else {
		die "Use either 'local', 'coverage' or 'global' for alignment type";
	}
	my @t = glob("$faDir/*.fa");
	my $faFiles = join( ' ', @t );
	print "$faFiles\n";

	my $cmd = "$aligner $faFiles -d $delim -o $alignDir/align.out.unsorted";
	print "Generating the identity scores ...\n";

	if ( system($cmd) ) {
		die "Cannot run $cmd\n";
	}
	if ($coverage == 1) {
		qx(sort < $alignDir/align.out.unsorted -n -k 5 -k 6 -t $delim | cut -d "$delim" -f 1,2,3,5 > $alignDir/align.out);
	} else {
		qx(sort < $alignDir/align.out.unsorted -n -k 5 -k 6 -t $delim | cut -d "$delim" -f 1,2,3,4 > $alignDir/align.out);
	}
	for ( my $hh = 1 ; $hh <= 12 ; $hh++ ) {
		my $cmd1 =
		  "$histogramer $faFiles -d $delim -k $hh -o $histDir/hist$hh.out";

		print "Generating the $hh-mers histogram ...\n";

		if ( system($cmd1) ) {
			die "Cannot run $cmd1\n";
		}
	}
}

sub runMatlab {
	my $cmd = "matlab -nodisplay -nosplash -r \"cd $matlabDir; EVALUATE(\'$outDir\'); exit;\"";
	print "Running MATLAB ...\n";
	print "sh -c $cmd\n";
	if ( system($cmd) ) {
		die "Could not run $cmd\n";
	}

}

sub getFreshDir {
	my ($dir) = @_;
	if ( -d $dir ) {
		print "Removing $dir ...\n";
		rmtree($dir) or die "Cannot remove $dir\n";
	}
	print "Making $dir ...\n";
	mkdir($dir) or die "Cannot make $dir\n";
}
