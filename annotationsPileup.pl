#!/usr/bin/perl

my $usage="
Annotations pileup: takes fasta file with contigs as main identifiers
(the first word after \'>\' symbol),
adds to the headers whatever it finds in the supplied space- or tab-delimited 
lookup tables:
contig2isogroup (c2i)
isogroup2genename (i2gn)
isogroup2GO (i2go)
isogroup2KEGG (i2k)

arguments:

required:
fasta=[file name]
c2i=[file name]

optional:
i2gn=[file name]
i2go=[file name]
i2k=[file name]

";

my $fasta;
if ("@ARGV"=~/fasta=(\S+)/) { $fasta=$1;}
else { die $usage;}
my $c2i;
if ("@ARGV"=~/c2i=(\S+)/) { $c2i=$1;}
else { die $usage;}
my $i2gn;
if ("@ARGV"=~/i2gn=(\S+)/) { $i2gn=$1;}
my $i2go;
if ("@ARGV"=~/i2go=(\S+)/) { $i2go=$1;}
my $i2k;
if ("@ARGV"=~/i2k=(\S+)/) { $i2k=$1;}

my $e1;
my $e2;
my @ee;

my %c2iso={};
open TT, $c2i or die "cannot open contigs to isogroups table $c2i\n";
while (<TT>) {
	chomp;
	($e1,$e2)=split(/\s/,$_);
	$c2iso{$e1}=$e2;
#print "$e1|$c2iso{$e2}|c2iso\n";
}
close TT;
my %i2gene={};
if ($i2gn){
	open TT, $i2gn or die "cannot open isogroups to genenames table $i2gn\n";
	while (<TT>) {
		chomp;
		($e1,@ee)=split(/\s/,$_);
#		$e2=join("_",@ee);
		$i2gene{$e1}="@ee";
#print "$e1|$i2gene{$e1}|i2gene\n";
	}
}
close TT;
my %i2GO={};
if ($i2go){
	open TT, $i2go or die "cannot open isogroups to GO table $i2go\n";
	while (<TT>) {
		chomp;
		($e1,$e2)=split(/\s/,$_);
		$i2GO{$e1}=$e2;
#print "$e1|$i2GO{$e1}|i2GO\n";
	}
}
close TT;
my %i2kegg={};
if ($i2k){
	open TT, $i2k or die "cannot open isogroups to KEGG table $i2k\n";
	while (<TT>) {
		chomp;
		($e1,$e2)=split(/\s/,$_);
		$i2kegg{$e1}=$e2;
#print "$e1|$i2kegg{$e1}|i2kegg\n";
	}
}
close TT;

open FF, $fasta or die "cannot open fasta file $fasta\n";
my $cname="";
my $seq="";
my @noiso=();
my $newname;

while (<FF>) {
	chomp;
	if ($_=~/^>(\S+)/) {
		$newname=$1;
		if ($seq && $cname) {
			my $iso="";
			if ($c2iso{$cname}) {
				$iso=$c2iso{$cname};
				print ">$cname gene=$iso";
			}
			else { 
				print ">$cname gene=NA ";
				push @noiso, $cname;
			}
			if ($i2gn && $i2gene{$iso}) {
				print " Gene=$i2gene{$iso}" unless ($i2gene{$iso}=~/annot missing/);	
			}
			if ($i2go && $i2GO{$iso}) {
				print " GOterms=$i2GO{$iso}";	
			}
			if ($i2k && $i2kegg{$iso}) {
				print " KEGGterms=$i2kegg{$iso}";	
			}
			print "\n$seq\n";
		}
		$cname=$newname;
		$seq="";
	}
	else {$seq.=$_;}
}
if ($seq && $cname) {
	my $iso="";
	if ($c2iso{$cname}) {
		$iso=$c2iso{$cname};
		print ">$cname gene=$iso";
	}
	else { 
		print ">$cname gene=NA ";
		$noiso++;
	}
	if ($i2gn && $iso && $i2gene{$iso}) {
		print " Gene=$i2gene{$iso}";	
	}
	if ($i2go && $iso && $i2GO{$iso}) {
		print " GOterms=$i2GO{$iso}";	
	}
	if ($i2k && $iso && $i2kegg{$iso}) {
		print " KEGGterms=$i2kegg{$iso}";	
	}
	print "\n$seq\n";
}

if ($noiso[0]) {
	print STDERR "\n",$#noiso+1," contigs without isogroup designation:
	@noiso 
	check whether your $c2i table is up to date
	";
}
	