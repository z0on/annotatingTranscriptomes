#!/usr/bin/perl

my $usage ="

Reformats the per-contig KOG class annotations into per-isogroup.
Removes isogroups whose contigs match different KOG classes.

Input:

   kogs=[filename] : output of getKOGs.pl
seq2iso=[filename] : contig to isogroup table

Prints to STDOUT

";

my $inkog;
my $s2i;
if ("@ARGV"=~/kogs=(\S+)/) { $inkog=$1;} else { die $usage;}
if ("@ARGV"=~/seq2iso=(\S+)/) { $s2i=$1;} else { die $usage;}

open INK, $inkog or die "cannot open $inkog\n";
open S2I, $s2i or die "cannot open $s2i\n";

my %s2iso={};
while (<S2I>) {
	chop;
	(my $s,my $i)=split("\t",$_);
	$s2iso{$s}=$i;
}

my %seen={};
my %len={};
my $i2k={};

while (<INK>) {
