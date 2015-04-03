#!/usr/bin/env perl
use lib "/home/01207/em23799/bin/BioPerl-1.6.0";
use Bio::Perl;
use Bio::SeqIO;

if ($ARGV[0] eq "-h") {print "usage: script input.fasta\n"; exit;}
unless ($#ARGV == 0) {print "usage: script input.fasta\n"; exit;}

my $infile = $ARGV[0];
my $seqs = new Bio::SeqIO(-file=>$infile, -format=>"fasta");

print"\n";
print $infile, "\n", "-"x25, "\n";
$countseq = 0;
$totbp = 0;
$ambigbp = 0;
$maxlen = 0;
$minlen = 1000000000000;
$countn = 0;
my %slh;
while($seq = $seqs->next_seq)
	{
	$lenbp = 0;
	$countseq++;
	$ss = $seq->seq;
	while($ss =~ /\w/g) {$lenbp++;}
	$slh{$seq->display_id} = $lenbp;
	foreach(split("",$ss)) 
		{
		if($_ !~ /[ACGT]/i){$ambigbp++;}
		if($_ =~ /N/i){$countn++;}
		}
	$totbp += $lenbp;		
	if($lenbp>$maxlen){$maxlen = $lenbp;}
	if($lenbp<$minlen){$minlen = $lenbp;}
	}
$meanbp = int($totbp/$countseq+0.5);
@ssia = sort{$slh{$a}<=>$slh{$b}}(keys(%slh));
my $tmpsum = 0;
foreach ($i=@ssia; $i>=0; $i--)
	{
	$tmpsum += $slh{$ssia[$i]};
	if ($tmpsum+$slh{$ssia[$i-1]}>=$totbp*0.5) {$n50=$slh{$ssia[$i]}; last;}
	}

print $countseq, " sequences.\n";
print $meanbp, " average length.\n";
print $maxlen, " maximum length.\n";
print $minlen, " minimum length.\n";
print "N50 = ", $n50, "\n";

if($totbp<1000)
	{
	print $totbp, " bp altogether.\n";
	print $ambigbp, " ambiguous bp. ";
	print "(", int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print $countn, " of those are Ns. ";
	print "(", int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
elsif ($totbp<1000000)	
	{
	print int($totbp/100+0.5)/10, " kb altogether ($totbp bp).\n";
	print int($ambigbp/100+0.5)/10, " ambiguous kb. ";
	print "(", $ambigbp, " bp, ";
	print int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print int($countn/100+0.5)/10, " kb of Ns. ";
	print "(", $countn, " bp, ";
	print int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
elsif ($totbp>=1000000)
	{
	print int($totbp/100000+0.5)/10, " Mb altogether ($totbp bp).\n";
	print int($ambigbp/100000+0.5)/10, " ambiguous Mb. ";
	print "(", $ambigbp, " bp, ";
	print int($ambigbp/$totbp*1000+0.5)/10, "%)\n";
	print int($countn/100000+0.5)/10, " Mb of Ns. ";
	print "(", $countn, " bp, ";
	print int($countn/$totbp*1000+0.5)/10, "%)\n";
	}
print "-"x25, "\n\n";

