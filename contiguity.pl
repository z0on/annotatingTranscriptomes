#!/usr/bin/perl

my $usage="

Calculates transcriptome contiguity as per
Martin and Wang, Nat Rev Genetics 12, 671-682, 2012:

Contiguity is the percentage of transcriptome-matching reference 
sequences aligned over [threshold] fraction of their length against 
the longest matching contig in the assembly.

For example, contiguity of 0.5 at the threshold 0.75 implies that
50% of reference sequences that were matched by the transcriptome
contigs were aligned over 75% or more of their length against
the longest matching contig.

Arguments:

        hits=[filename] : *_hits.tab output file from CDS_extractor_v2.pl 
threshold=[from 0 to 1] : aligned length threshold for contiguity calculation; 
                          default 0.75
                          
Output:

prints the calculated contiguity to STDOUT
generates a pdf historgam of aligned proportions (calls Rscript)

Mikhail Matz, matz\@utexas.edu
October 2014

";

my $hits;
if ("@ARGV"=~/hits=(\S+)/) { $hits=$1; } else { die $usage;}
my $thresh=0.75;
if ("@ARGV"=~/threshold=(\S+)/) { $thresh=$1; } 

open HITS, $hits or die "cannot open hits file $hits\n";

my %aligned={};
my @hline;
my $fl=1;

while (<HITS>) {
	if ($fl==1) { $fl-- and next;}
	chomp;
	@hline=split(/\t/,$_);
	my $cov=sprintf("%.3f",$hline[2]/$hline[4]);
	my $hit=$hline[3];
	push @{$aligned{$hit}}, $cov;
}

my @hnames=keys(%aligned);

my %bestalig={};
my $contig;

my $rfname=$hits."_contiguity.R";
open RF, ">$rfname";

print {RF} "aligns=c(";

foreach $hit (@hnames) {
	next if ($hit=~/HASH/);
	my $maxcov=0;
	foreach $cov (@{$aligned{$hit}}) {
		if ($cov>$maxcov) { 
			if ($cov>1.3){$cov=1.3;}
			$maxcov=$cov;
			$bestalig{$hit}=$cov;
		}
	}
	if ($bestalig{$hit}>=$thresh) { $contig++;}
	if ($hit eq $hnames[0]) { print {RF} "$bestalig{$hit}";}
	else { print {RF} ",$bestalig{$hit}";}
}

$contig=sprintf("%.2f",$contig/($#hnames+1));
print "\ncontiguity at $thresh threshold: $contig\n\n";
print {RF} ")

pdf(\"contiguity_$hits.pdf\", width=4,height=4.5)
hist(aligns,xlim=c(-0.2,1.3),breaks=30,
	bty=\"n\",
	main=paste(\"Contiguity at \",$thresh,\" threshold: \",$contig,sep=\"\"),
	xlab=\"reference CDS coverage\",
	mgp=c(2.3,1,0)
)
dev.off()
";
close RF;
system("Rscript $rfname");
