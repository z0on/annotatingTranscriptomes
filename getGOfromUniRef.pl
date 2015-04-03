#!/usr/bin/perl

my $usage="

Assigns annotations to transcriptome (per-isogroup), based on best blast
hit per isogroup. Fasta headers of the transcriptome may or may not contain 
'gene=isogroupNNN' label, if not, the annotation will be made for the individual
sequence.

Extracts GO terms from uniProt's idmapping_selected.tab file, 
which must be present in the same directory
(wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz)

assumes that we are workign with uniref50

Outputs two tab-delimited tables: iso2go (GO annotations) and iso2hit (hit identifiers)

Arguments:

blast=<file name> name of the blast file
prefix=<text> prefix for the output file names
fastaQuery=<file name> fasta file of the transcriptome, to extract isogroups
(optional) evalue=<float> evalue cutoff, default 0.001


";

if ("@ARGV"=~/blast=(\S+)/) { $inp=$1;} else { die $usage;}
if ("@ARGV"=~/prefix=(\S+)/) { $pre=$1;} else { die $usage;}
if ("@ARGV"=~/fastaQuery=(\S+)/) { $fas=$1;} else { die $usage;}
$fasdb="idmapping_selected.tab";
my $ecut=0.001;
if ("@ARGV"=~/evalue=(\S+)/) { $ecut=$1;} 

my %iso2gene={};
my %iso2hit={};
my $iso;
my $hit="";
my $desc="";
my $ev="";
my $add=0;

warn "reading fasta..\n";
open FAS, $fas or die "cannot open fasta file $fas\n";
my %f2iso={};
while (<FAS>){
	chop;
	if ($_=~/^>(\S+).+gene=(isogroup\d+)/){ $f2iso{$1}=$2;}
	elsif ($_=~/^>(\S+)/){ $f2iso{$1}=$1;}
#print "$1\t$f2iso{$1}\n";
}
close FAS;

warn "reading database...\n";
open FAS, $fasdb or die "cannot open uniRef idmapping file $fasdb\n";
my %fdb2gname={};
while (<FAS>){
	chop;
	my $gos;
	if ($_=~/\t(GO.+)\tUniRef100/){ 		
		$gos=$1;
		$gos=~s/\s//g;
	}
	if ($_=~/(UniRef50\S+)/){ 
		$fdb2gname{$1}=$gos;
print "entry: $1  GO:  $gos\n";
	}
}
close FAS;

warn "processing blast file...\n";
open INP, $inp or die "cannot open blast file $inp\n";
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
		if ($line=~/undetermined/) { next;}
#print "eline: $line\n$1\n";
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
close INP;

#$outhits=$pre."_iso2hit.tab";
$outgenes=$pre."_iso2go.tab";
#open HI, ">$outhits" or die "cannot create $outhits\n";
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
#	print {HI} "$iso\t${$iso2hit{$iso}}[$eind]\t${$evals{$iso}}[$eind]\n";
	print {GE} "$iso\t${$iso2gene{$iso}}[$eind]\n";
}
	
	
	
	
	
	
	
	
	
	