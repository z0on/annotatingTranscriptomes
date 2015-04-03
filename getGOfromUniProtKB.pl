#!/usr/bin/perl

my $usage="

Assigns annotations to transcriptome (per-isogroup), based on best blast
hit per isogroup. Fasta headers of the transcriptome may or may not contain 
'gene=isogroupNNN' label, if not, the annotation will be made for the individual
sequence.

Extracts GO terms from uniProt's idmapping_selected.tab file, 
which must be present in the same directory
(wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz)

assumes that the blast was to uniprot_sprot.fasta (uniprot swssprot KB)

Outputs a tab-delimited table <prefix>_iso2go.tab 

Arguments:

blast=<file name> name of the blast file
prefix=<text> prefix for the output file names
fastaQuery=<file name> fasta file of the transcriptome, to extract isogroups
(optional) evalue=<float> evalue cutoff, default 0.0001

";

if ("@ARGV"=~/blast=(\S+)/) { $inp=$1;} else { die $usage;}
if ("@ARGV"=~/prefix=(\S+)/) { $pre=$1;} else { die $usage;}
if ("@ARGV"=~/fastaQuery=(\S+)/) { $fas=$1;} else { die $usage;}
$fasdb="idmapping_selected.tab";
my $ecut=0.0001;
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

my %usedDesc={};
warn "processing blast file...\n";
open INP, $inp or die "cannot open blast file $inp\n";
while(<INP>){
	chop;
	my $line=$_;
	if ($line=~/Query=\s+(\S+)/) {
		if ($hit && $ev<=$ecut) { 
			push @{$iso2hit{$iso}},$hit;
			push @{$evals{$iso}},$ev;
		}
		$iso=$f2iso{$1};
#print "----------------\n$iso\n";
		$hit="";
		$ev="";
	}
	elsif ($line=~/^\s+(\S+).+\s(\d\S+\d)$/ && !$hit) {
		if ($line=~/undetermined/) { next;}
		my @hitparts=split(/\|/,$1);
		$hit=$hitparts[1];
		$usedDesc{$hit}=1;
		$ev=$2;
#print "----------------\n$iso\n";
#print "hit:$hit E:$ev Gene:$desc\n";
	}
}
if ($hit && $ev<=$ecut) { 
	push @{$iso2hit{$iso}},$hit;
	push @{$evals{$iso}},$ev;
}
close INP;


warn "reading database...\n";
open FAS, $fasdb or die "cannot open uniprot idmapping file $fasdb\n";
my %hit2go={};
while (<FAS>){
	chop;
	my $gos;
	my $ggn;
	if ($_=~/^(\S+)\s.+\t(GO.+?)\t/){
		$ggn=$1;
		$gos=$2;
		next unless ($usedDesc{$ggn});
		$gos=~s/\s//g;
		$hit2go{$ggn}=$gos;		
#print "entry:$ggn|GO::$hit2go{$ggn}\n";
	}
}
close FAS;

$outgenes=$pre."_iso2go.tab";
open GE, ">$outgenes" or die "cannot create $outgenes\n";

foreach $iso (keys %iso2hit) {
	next if ($iso=~/HASH/);
	my $emin=10;
	my $eind=0;
#print "\t\t${$iso2hit{$iso}}[$e]\t${$evals{$iso}}[$e]\n";
	if (${$evals{$iso}}[1]) {
		for ($e=0;$ev=${$evals{$iso}}[$e];$e++){
			if ($ev<$emin) {
				$eind=$e;
				$emin=$ev;
			}
		}
	}
	$hit=${$iso2hit{$iso}}[$eind];
	print {GE} "$iso\t$hit2go{$hit}\n";
}
