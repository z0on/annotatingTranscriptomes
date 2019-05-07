#!/usr/bin/perl

$usage= "

fasta2SBH: prepares a transcriptome assembly for annotation using single-best-hit 
(SBH) method of the KEGG KAAS tool, by throwing out all isoforms except the longest one 
for each gene.

arg1: transcriptome fasta file. Assumes there is 'isogroup' annotations in the
	  fasta header for 454 assemblies, such as 'isogroup00234', or the fasta header has 
	  form  'c36966_g1_i1' (Trinity). If not, will use the first 
	  word in the fasta header as the group identifier (which will not result in any 
	  slimming-down since those are supposed to be unique).
	  
output: prints the slimmed-down fasta data to sdout.

";

open db, $ARGV[0] or die $usage; 

my %ig2seq={};
my %ig2def={};
my $def;
my $seq;
my $ig;
my $tig;
while (<db>){
	$line=$_;
	if ($line=~/^>/){
		if ($seq){
			if (length($ig2seq{$ig})<length($seq)) { 
				$ig2seq{$ig}=$seq;
				$ig2def{$ig}=$def;
			}
		}
		if ($line=~/(isogroup\S+)/) { $ig=$1;}
		elsif ($line=~/^>(\S+\d+_g\d+)_i\d+/) { $ig=$1; }
		elsif ($line=~/^>(\S+)/) { $ig=$1; }
		else { die "stopped: can't decipher fasta header: \n$line\n"; } 
		$def=$line;
		$seq="";
	}
	else { 
		chomp($line);
		$seq.=$line;
	}
}	

foreach $ig (keys %ig2seq){
	next if ($ig=~/HASH/);
	if (!$ig2def{$ig}) { die "sequence $ig does not have a header!\n";}
	my $header=substr $ig2def{$ig},1;
	print ">$ig ",$header,$ig2seq{$ig},"\n";
}
