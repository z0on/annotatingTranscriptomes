#!/usr/bin/perl

my $usage="

Assigns annotations to transcriptome (per-isogroup), 
based on match to KOG database generated at WebMGA:
http://weizhong-lab.ucsd.edu/metagenomic-analysis/server/kog/ 

Fasta headers of the transcriptome may or may not contain 
'gene=GeneID' label (for Trinity, GeneID is the component), 
if not, the annotation will be made for the individual
contigs.

Outputs two tab-delimited tables: 
iso2kogID
iso2kogDef (KOG deflines - like gene names)
iso2kogClass (KOG classes)

Arguments:

fastaQuery=<file name> fasta file of the transcriptome, to extract isogroups
prefix=<text> prefix to put into output file names
kogMatch=<file name> KOG search result (output.2 from WebMGA))
(optional) evalue=<float> evalue cutoff, default 0.0001

";

if ("@ARGV"=~/prefix=(\S+)/) { $pre=$1;} else { die $usage;}
if ("@ARGV"=~/fastaQuery=(\S+)/) { $fas=$1;} else { die $usage;}
if ("@ARGV"=~/kogMatch=(\S+)/) { $fasdb=$1;} else { die $usage;}
my $ecut=0.0001;
if ("@ARGV"=~/evalue=(\S+)/) { $ecut=$1;} 

my %iso2gene={};
my %iso2kogID={};
my %iso2kogClass={};
my %iso2hit={};
my $iso;
my $hit="";
my $desc="";
my $ev="";
my $add=0;

open FAS, $fas or die "cannot open fasta file $fas\n";
my %f2iso={};
while (<FAS>){
	chop;
	if ($_=~/^>(\S+).+gene=(\S+)/){ $f2iso{$1}=$2;}
	elsif ($_=~/^>(\S+)/){ $f2iso{$1}=$1;}
	
#print "$1\t$f2iso{$1}\n";
}


open FAS, $fasdb or die "cannot open KOG match file $fasdb\n";
my %fdb2gname={};
my %fdb2kogID={};
my %fdb2kogClass={};
my $kc="";
while (<FAS>){
	chop;
	@ld=split(/\t/,$_);
	if ($ld[1]!~/KOG/) { next;}
	$kc=$ld[$#ld];
	$kc=~s/ $//;
	if ($kc!~/General function prediction only|Function unknown/ && $ld[2]<$ecut && "@{$fdb2kogClass{$f2iso{$ld[0]}}}"!~/$kc/){
		push @{$fdb2kogClass{$f2iso{$ld[0]}}}, "$kc" ;
		push @{$fdb2kogID{$f2iso{$ld[0]}}}, $ld[1];
		push @{$fdb2gname{$f2iso{$ld[0]}}}, "$ld[$#ld-2]";
	}
}

$outgenes=$pre."_iso2kogDef.tab";
$outkog=$pre."_iso2kog.tab";
$outkogClass=$pre."_iso2kogClass.tab";
open GE, ">$outgenes" or die "cannot create $outgenes\n";
open KOG, ">$outkog" or die "cannot create $outkog\n";
open KC, ">$outkogClass" or die "cannot create $outkogClass\n";

foreach $iso (keys %fdb2gname) {
	next if ($iso=~/HASH/);
	for (my $kk=0;${$fdb2kogID{$iso}}[$kk];$kk++) {
		print {KC} "$iso\t${$fdb2kogClass{$iso}}[$kk]\n";
		print {GE} "$iso\t${$fdb2gname{$iso}}[$kk]\n";
		print {KOG} "$iso\t${$fdb2kogID{$iso}}[$kk]\n";
	}
}
