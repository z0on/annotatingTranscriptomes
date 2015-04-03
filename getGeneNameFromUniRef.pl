#!/usr/bin/perl

my $usage="

Assigns annotations to transcriptome (per-isogroup), based on best blast
hit per isogroup. Fasta headers of the transcriptome may or may not contain 
'gene=isogroupNNN' label, if not, the annotation will be made for the individual
sequence.

Extracts gene names from fasta headers of the hits.

Outputs two tab-delimited tables: iso2gene (gene names) and iso2hit (hit identifiers)

Arguments:

blast=<file name> name of the blast file
prefix=<text> prefix to put into output file names
fastaQuery=<file name> fasta file of the transcriptome, to extract isogroups
(optional) evalue=<float> evalue cutoff, default 0.001


";

if ("@ARGV"=~/blast=(\S+)/) { $inp=$1;} else { die $usage;}
if ("@ARGV"=~/prefix=(\S+)/) { $pre=$1;} else { die $usage;}
if ("@ARGV"=~/fastaQuery=(\S+)/) { $fas=$1;} else { die $usage;}
$fasdb="uniref50.fasta";
my $ecut=0.001;
if ("@ARGV"=~/evalue=(\S+)/) { $ecut=$1;} 

open INP, $inp or die "cannot open blast file $inp\n";

my %iso2gene={};
my %iso2hit={};
my $iso;
my $hit="";
my $desc="";
my $ev="";
my $add=0;

warn "reading query fasta...\n";
open FAS, $fas or die "cannot open fasta file $fas\n";
my %f2iso={};
while (<FAS>){
	chop;
	if ($_=~/^>(\S+).+gene=(isogroup\d+)/){ $f2iso{$1}=$2;}
	elsif ($_=~/^>(\S+)/){ $f2iso{$1}=$1;}
# print "$1\t$f2iso{$1}\n";
}

warn "reading uniRef50 fasta...\n";
open FAS, $fasdb or die "cannot open fasta file $fasdb\n";
my %fdb2gname={};
while (<FAS>){
	chop;
	if ($_=~/^>(\S+)\s+(.+)\sn=/){ 
		$fdb2gname{$1}=$2 unless ($2=~/ncharacterized|ypothetical/);
#print "$1\t$fdb2gname{$1}\n";
	}
}
close FAS;

warn "processing blast...\n";
while(<INP>){
	chop;
	my $line=$_;
	if ($line=~/Query=\s+(\S+)/) {
		if ($hit && $ev<$ecut && $desc) { 
			push @{$iso2gene{$iso}},$desc;
			push @{$iso2hit{$iso}},$hit;
			push @{$evals{$iso}},$ev;
		}
		$iso=$f2iso{$1};
#print "----------------\n$iso\n";
		$hit="";
		$desc="";
		$ev="";
	}
	elsif ($line=~/^\s+(\S+).+\s(\d\S+\d)$/ && !$hit) {
		if ($line=~/ndetermined|ncharacterized|ypothetical|redicted protein/) { next;}
#print "eline: $line\n$1\t$2\n";
		$hit=$1;
		$ev=$2;
		$desc=$fdb2gname{$1};
#if ($hit=~/iso/) { 
#print "----------------\n$iso\n";
#print "hit:$hit E:$ev Gene:$desc\n";
#}
	}
}
if ($hit && $ev<=$ecut && $desc) { 
	push @{$iso2gene{$iso}},$desc;
	push @{$iso2hit{$iso}},$hit;
	push @{$evals{$iso}},$ev;
}

$outhits=$pre."_iso2hit.tab";
$outgenes=$pre."_iso2gene.tab";
open HI, ">$outhits" or die "cannot create $outhits\n";
open GE, ">$outgenes" or die "cannot create $outgenes\n";

foreach $iso (keys %iso2gene) {
	next if ($iso=~/HASH/);
	my $emin=10;
	my $eind=0;
	if (${$evals{$iso}}[1]) {
		for ($e=0;$ev=${$evals{$iso}}[$e];$e++){
# print "\t\t${$iso2hit{$iso}}[$e]\t${$evals{$iso}}[$e]\n";
			if ($ev<$emin) {
				$eind=$e;
				$emin=$ev;
			}
		}
	}
	print {HI} "$iso\t${$iso2hit{$iso}}[$eind]\t${$evals{$iso}}[$eind]\n";
	print {GE} "$iso\t${$iso2gene{$iso}}[$eind] E(blastx)=${$evals{$iso}}[$eind]\n";
}
	
	
	
	
	
	
	
	
	
	