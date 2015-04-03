#!/usr/bin/perl

my $usage= "
addAnnotField:

inserts \"gene=QWERTY\" annotations into fasta headers 
where QUERTY is the gene name (or cluster ID) taken from
a tab-delimited table of seq id - gene name

Arguments: 
1. fasta file
2. seqID-geneName table

prints to STDOUT.

";

my $fasta=shift(@ARGV) or die $usage;
my $table=shift(@ARGV) or die $usage;

open TAB, $table or die "$usage\n\ncannot open gene names table $table\n";
open FAS, $fasta or die "$usage\n\ncannot open fasta $fasta\n";

my %seq2gene={};
my $seq;
my $gene;
while (<TAB>) {
	($seq,$gene)=split(/[\t,]/,$_);
	$seq=~s/\.([\w\d]).+/\.$1/;
	$seq2gene{$seq}=$gene;
}
close TAB;

while (<FAS>) {
	chop;
	if ($_=~/^>(\S+)/) {
		$seq=$1;
		$seq=~s/\.([\w\d]).+/\.$1/;
		if ($seq2gene{$seq}) { print "$_ gene=$seq2gene{$seq}\n";}
		else { 
warn "no gene for $seq\n";
			print "$_ gene=$seq\n";
			}
	}
	else { print "$_\n";}
}
close FAS;