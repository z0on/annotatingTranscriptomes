#!/usr/bin/perl

my $usage="
restores pairing among quality-filtered paired-end fastq files
Writes two fully paired files and one with all the unpaired reads.

arg1: left reads
arg2: right reads

output: 
R1_[left reads filename] - right paired reads
R2_[right reads filename] - left paired reads
Unp_[left reads filename]_[right reads filename] - unpaired reads

";

my $lf=shift or die $usage;
my $rf=shift or die $usage;

my %rights={};
my %paired=();

open RR,$rf or die "cannot open right reads $rf\n";
my $count=0;
while (<RR>) {
	$count++;
	if ($count==5) { $count=1;}
	if ($count==1){
		chomp;
		$_=~s/(\S+)\s.+/$1/;
		$rights{$_}=1;
	}
}
close RR;

open LL,$lf or die "cannot open left reads $lf\n";
my $count=0;
while (<LL>) {
	$count++;
	if ($count==5) { $count=1;}
	if ($count==1){
		chomp;
		$_=~s/(\S+)\s.+/$1/;
		if ($rights{$_}==1) { 
			$paired{$_}=1;
		}
	}
}
undef %rights;
close LL;


my $unname="Unp_".$lf."_".$rf;
my $lname="R1_".$lf;
my $rname="R2_".$rf;
open UN,">$unname";

open RO,">$rname";
open RR,$rf;
my $count=0;
my $pri=0;
while (<RR>) {
	chomp;
	$count++;
	if ($count==5) { $count=1;}
	if ($count==1){
		$_=~s/(\S+)\s.+/$1/;
		if ($paired{$_}){
			$pri=1;
		}
		else { 
			$pri=0;
		}	
	}
	if ($pri==1) { print {RO} $_,"\n"; }
	else { print {UN} $_,"\n"; }
}
close RR;
close RO;

open LO,">$lname";
open LR,$lf;
my $count=0;
my $pri=0;
while (<LR>) {
	chomp;
	$count++;
	if ($count==5) { $count=1;}
	if ($count==1){
		$_=~s/(\S+)\s.+/$1/;
		if ($paired{$_}){
			$pri=1;
		}
		else { $pri=0;}	
	}
	if ($pri) { print {LO} $_,"\n"; }
	else { print {UN} $_,"\n"; }
}
close LR;
close LO;
close UN;