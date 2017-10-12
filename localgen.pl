#!/usr/bin/env perl
#
# Author: Benjamin T James
#

use warnings;
use strict;
use Bio::SeqIO;
use List::Util qw(shuffle);

my $num = 0;
my ($file) = @ARGV;
die "Usage: $0 input.fasta" unless $#ARGV == 0;

my $seqio = Bio::SeqIO->new(-file => $file, '-format' => 'Fasta');

my @db = ();
while (my $seq = $seqio->next_seq) {
	push @db, $seq;
	$num++;
}

my @indices = shuffle(0..$#db);
my @sample = @db[ @indices ];

print ">" . $sample[0]->display_id . "|template\n";
print $sample[0]->seq . "\n\n";

my $template = $sample[0]->seq;
my $half_num = int($num / 2);

for (my $i = 1; $i < $half_num; $i++) {
	my $first = lc($sample[$i]->seq);
	my $second = rand_substr($template);

	my $len = length $first;
	my $begin = rand_substr(substr($first, 0, int(2 * $len / 3)));
	my $end = rand_substr(substr($first, int($len / 3)));

	print ">" . $sample[$i]->display_id . "|$i|split\n";
	print "$begin$second$end\n\n";
	print ">" . $sample[$i]->display_id . "|" . $i . "_actual\n";
	print "$first\n\n";
}

for (my $i = $half_num; $i < $num; $i++) {
	my $first = lc($sample[$i]->seq);
	my $second = rand_substr($template);

	my $len = length $first;
	my $begin = rand_substr(substr($first, 0, int(2 * $len / 3)));
	my $end = rand_substr(substr($first, int($len / 3)));

	print ">" . $sample[$i]->display_id . "|$i|prefixed\n";
	print "$begin$second\n\n";
	print ">" . $sample[$i]->display_id . "|$i|suffixed\n";
	print "$second$end\n\n";
	print ">" . $sample[$i]->display_id . "|" . $i . "_actual|split\n";
	print "$first\n\n";
}

sub rand_substr {
	my $str = shift;
	my $len = length($str);
	my $newlen = int(rand($len));
	my $newoffset = int(rand($len - $newlen));
	return substr($str, $newoffset, $newlen);
}
