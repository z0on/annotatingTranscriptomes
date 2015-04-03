#!/usr/bin/perl
my $usage="

Removes duplicated reads prior to transcriptome assembly.
Input must be fastq files, trimmed, quality-filtered and 
sorted into left, right and unpaired reads.

Identifies duplicates based on identity of bases 5-30
on both left and right ends (if paired reads are supplied)
or in unpaired reads.

If all three types are given, unpaired are processed last
and are gettivn discarded if they find a match either 
among left or right reads

Usage:

dedupTranscriptome.pl left=[filename] right=[filename] unp=[filename]

output: same file names with .dedup extention

Mikhail Matz, August-November 2014
matz\@utexas.edu

";


my $left="";
my $right="";
my $unp="";

my $rlcount=0;
if ("@ARGV"=~/left=(\S+)/) { 
	$left=$1;
	$rlcount++;
	open LL,$left or die "cannot open $left\n";
	my $outleft=$left.".dedup";
	open LLO,">$outleft";
}
if ("@ARGV"=~/right=(\S+)/) { 
	$right=$1;
	$rlcount++;
	open RR,$right or die "cannot open $right\n";
	my $outright=$right.".dedup";
	open RRO,">$outright";
}
if ("@ARGV"=~/unp=(\S+)/) { 
	$unp=$1;
	open UU,$unp or die "cannot open $unp\n";
	my $outunp=$unp.".dedup";
	open UUO,">$outunp";
}
if ($rlcount >0 ){
	if ($rlcount<2){ die "$usage
--------------------------
specify either both left and right reads 
plus unpaired, or unpaired only\n";}
}
elsif (!$unp) {die "$usage
--------------------------
specify either both left and right reads 
plus unpaired, or unpaired only\n";}

my %seenleft={};
my %seenright={};
my %seenunp={};

my $name="";
my $seql="";
my $seqr="";
my $ql="";
my $qr="";
my $ii=1;
my $lrcount=0;
my $lrall==0;

if ($left && $right) {
	while (<LL>) {
		chop;
		my $llin=$_;
		$rlin=<RR>;
		chop($rlin);
		if ($ii==1 ) {
			$lrall++;
			if (!$qr || !$ql) {
				$name=$llin;
				$ii=2;
				next;
			}
			if (!$seenleft{substr($seql,4,25)} or !$seenright{substr($seqr,4,25)}) {
				$seenleft{substr($seql,4,25)}=1;
				$seenright{substr($seqr,4,25)}=1;
				$lrcount++;
				print {LLO} "$name\n$seql\n+\n$ql\n";
				print {RRO} "$name\n$seqr\n+\n$qr\n";
			}
			$name=$llin;
			$seql="";
			$seqr="";
			$ql="";
			$qr="";
			$ii=2;
		}
		elsif ($ii==2) {
			$seql=$llin;
			$seqr=$rlin;
			$ii=3;
		}
		elsif ($ii==3){ 
			$ii=4;
		}
		else {
			$qr=$rlin;
			$ql=$llin;
			$ii=1;
		}
	}
	if (!$seenleft{substr($seql,4,25)} or !$seenright{substr($seqr,4,25)}) {
		$seenleft{substr($seql,4,25)}=1;
		$seenright{substr($seqr,4,25)}=1;
		$lrcount++;
		print {LLO} "$name\n$seql\n+\n$ql\n";
		print {RRO} "$name\n$seqr\n+\n$qr\n";
	}
}
print "PAIRED: kept $lrcount out of $lrall\n";
close LL;
close RR;
close LLO;
close RRO;

my $sequ="";
my $qu="";
my $ii=1;
my $ucount=0;
my $uall==0;

if ($unp) {
	while (<UU>) {
		chop;
		my $ulin=$_;
		if ($ii==1 ) {
			$uall++;
			if (!$qu) {
				$name=$ulin;
				$ii=2;
				next;
			}
			if (!$seenleft{substr($sequ,4,25)} && !$seenright{substr($sequ,4,25)} && !$seenu{substr($sequ,4,25)}) {
				$seenu{substr($sequ,4,25)}=1;
				$ucount++;
				print {UUO} "$name\n$sequ\n+\n$qu\n";
			}
			$name=$ulin;
			$sequ="";
			$qu="";
			$ii=2;
		}
		elsif ($ii==2) {
			$sequ=$ulin;
			$ii=3;
		}
		elsif ($ii==3){ 
			$ii=4;
		}
		else {
			$qu=$ulin;
			$ii=1;
		}
	}
	if (!$seenleft{substr($sequ,4,25)} && !$seenright{substr($sequ,4,25)} && !$seenu{substr($sequ,4,25)}) {
		$seenu{substr($sequ,4,25)}=1;
		$ucount++;
		print {UUO} "$name\n$sequ\n+\n$qu\n";
	}
}
print "UNPAIRED: kept $ucount out of $uall\n";