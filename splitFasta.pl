#!/usr/bin/perl

my $fas=shift or die "specify fasta file to split\n";
my $nc = shift or die "specify number of chunks\n";

my $l=`grep ">" $fas | wc -l`;
chop $l;
my $chunk=sprintf("%.0f",$l/$nc);
if ($chunk*$nc<$l) { $chunk++;}
print "splitting $l records into $nc chunks; $chunk records per chunk\n"; 
my $k=1;
my $outname="subset".$k."_".$fas;
open INP, $fas or die "cannot open input $fas\n";
open OUT, ">$outname" or die "cannot create $outname\n";
my $counter=0;
while (<INP>) {
	if ($_=~/^>/) { 
		$counter++;
		if ($counter>$chunk) {
			$k++;
			print "chunk",$k-1," ",$counter-1," seqs\n";
			$counter=1;
			$outname="subset".$k."_".$fas;
			open OUT, ">$outname" or die "cannot create $outname\n";
		}
	}
	print {OUT} $_;
}
print "chunk",$k," ",$counter," seqs\n";

