#!/usr/bin/perl

$ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
$ENV{BIOPERL_INDEX} = ".";
use Bio::DB::Fasta;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO; 
#use Bio::DB::GenBank;

if (!$ARGV[1]) {
print "

 CDS_extractor_v2.pl (version 2.2 , September 10 2014)

 extracts CDS (protein coding sequences) from a messy next-gen assembly and translates 
 them based on blastx results; corrects frameshifts due to indels

 Requires Bioperl ( http://www.bioperl.org/wiki/ )

 Usage: CDS_extractor_v2.pl [fasta file] [blastx results file] 
 optional additional arguments :
 	allhits    : will cause the script to consider matches to predicted AND 
 	             hypothetical proteins
	bridgegaps : will cause the script to translate linkers between neighboring 
	             in-frame HSPs, if there are no stops
	stopcheck  : (on by default,specify stopcheck=no to switch off) check extracted 
	             ORFs for stop codons, clip ends if a stop found near ORF end
	verbose    : print CDS extraction info to STDOUT
 
 Output:
 	[fasta filename]_PRO.fasta : translations based strictly on merged HSPs
 	[fasta filename]_CDS.fasta : corresponding coding sequences
 	[fasta filename]_PROends.fasta : have extended HSP ends until the nearest stop codon
 	[fasta filename]_CDSends.fasta : corresponding coding sequence
 	[fasta filename]_hits.tab : table of hits and their lengths
 
 Mikhail V Matz 
 University of Texas at Austin
 matz\@utexas.edu
 
\n\n";
exit(0);
}
#default parameters:
my $hsp_eval=0.0001 ; # e-vaue cutoff for HSPs
my $allhits=0; # use matches to proteins that are hypothetical, predicted etc.
my $stopcheck=1; # check extracted ORFs for stop codons, clip ends if one found near end
my $ends=1; # produce additional files with ORFs extended both ways to the nearest stop codon
my $linktrans=0; # translate the gap between neighboring HSPs if they are in the same frame and have no stops in between 
my $verbose=0;

$argline=join(" ",@ARGV);
if ($argline=~/ allhits/) { $allhits=1;}
if ($argline=~/ brigdegaps/) { $linktrans=1;}
if ($argline=~/ verbose/) { $verbose=1;}
if ($argline=~/ stopcheck=no/) { $stopcheck=0;}

if (!$ARGV[1]) { die "Usage: CDS_extractor.pl [fasta file to translate] [blastx results file] \n";}
my $prefix=$ARGV[0];
$prefix=~s/(\S+)\..+$/$1/;
my $proname=$prefix."_PRO.fas";
my $cdsname=$prefix."_CDS.fas";
my $proname1=$prefix."_PROends.fas";
my $cdsname1=$prefix."_CDSends.fas";

#-------------
#creating index, opening output files 

my $inx = Bio::DB::Fasta->new($ARGV[0]);

$outdna = Bio::SeqIO->new(-file => ">$cdsname",
                       -format => 'Fasta'
                        -alphabet => 'dna');
$outpro = Bio::SeqIO->new(-file => ">$proname",
                       -format => 'Fasta'
                       -alphabet => 'protein');
$outdnaE = Bio::SeqIO->new(-file => ">$cdsname1",
                       -format => 'Fasta'
                        -alphabet => 'dna');
$outproE = Bio::SeqIO->new(-file => ">$proname1",
                       -format => 'Fasta'
                       -alphabet => 'protein');
                       

#$in= Bio::SeqIO->new(-file => $ARGV[0],
#                           -format => 'Fasta');

my $tabname=$prefix."_hits.tab";
open HITS, ">$tabname" or die "cannot create $tabname\n";
print {HITS} "query\tqlen.bp\tqtrans.aa\thit\thlen.aa\n";

#-----------
# reading BLAST output

use Bio::SearchIO; 
my $inb = new Bio::SearchIO(-format => 'blast',  -file   => $ARGV[1]);
while( my $result = $inb->next_result ) {
	my $seqobj=$inx->get_Seq_by_id($result->query_name);
if ($verbose) {    print "--------------\nQuery=",   $seqobj->display_id,"\n"; }
    my $qlength=$result->query_length;
#	print $seqobj->seq(), "\n";

#---------
# looking for the best hit:  well-annotated, making the evalue cutoff, and
# showing the highest number of conserved residues in non-overlapping portions of HSPs

	my $closest="";
	my %coordinates={};
	my $besthit="";
	my $reverse=0;
	my %coord={};
	my %hsps={};
	my $bh_length=0;
		
  	while( my $hit = $result->next_hit ) {
  		$hitname=$hit->name;
    	if ($hit->description=~/predicted|hypothetical|similar|novel|unnamed/i){
#    		print "skipping hit:\n", $hit->description, "\n $1 \n";
    		next unless ($argline=~/ allhits/);
    	}
#    	next if ($hit->significance>$ARGV[2]);
    	next if ($hit->num_hsps=="-");
		my $covered=0;
		my $conserved=0;
		my $rcom=0;
		my @ends=();
		my %start={};
		my $hlength=$hit->length;
#print "Hit: $hitname\n";

		while( my $hsp = $hit->next_hsp ) {
#print "\tHSP: ", $hsp->start('query') . "_". $hsp->end('query'), "\n";
			if ($hsp->significance<$hsp_eval){
				$newend=$hsp->end('query');
				$newstart=$hsp->start('query');
				push @ends, $newend;
				$start{$newend}=$newstart unless (exists($start{$newend}));
				push @{$coord{$hitname}}, $hsp->start('query'), $hsp->end('query');
				push @{$hsps{$hitname}}, $hsp->start('query') . "_" . $hsp->end('query');
				if ($hsp->query->strand < 0) {$rcom=1;}
				else {$rcom=0;}
				$conserved+=$hsp->frac_conserved;
#print "\t",$hsp->start('query'), " - ", $hsp->end('query'), "\tconserv: ", $hsp->frac_conserved, "\tstrand: ", $hsp->query->strand, "\n";
			}
    	}
    	
    	$conserved=$conserved/$hit->num_hsps; # calculating average fraction of conserved residues
#print "\taverage conserved: $conserved\n";
    	my @coordinates=sort {$a <=> $b} @{$coord{$hitname}};
		for (my $h=0 ; $coordinates[$h] ; $h=$h+2) { # calculating the length of query covered without overlaps
			$covered+=$coordinates[$h+1]-$coordinates[$h];
		}
		$similar=$covered*$conserved; # key test statistic: approximate number of conserved residues in non-overlapping hsps
    	if ($similar>$closest){ 
    		$closest=$similar;
    		$besthit=$hitname;
    		$reverse=$rcom;
    		$bh_length=$hlength;
    	}
    }
if ($verbose) {    	print "best hit: $besthit  revcom:$reverse\n";  }

# -------------------
# merging overlapping HSPs that are in the same frame

	my @rehsp = ();
	my @recoord=();
	my @merge = ();
	@{$hsps{$besthit}}=sort {$a <=> $b} @{$hsps{$besthit}};
if ($verbose) {    print "raw HSPs: ", join(" ",@{$hsps{$besthit}}), "\n"; }
	for ($h=0; $hspp1=${$hsps{$besthit}}[$h];$h++) {
		if ($merge[$h]) {next;}
		($s1,$e1)=split(/_/,$hspp1);
#print "$hspp1 merge:\n";
		my $merged=0;
		for ($h2=$h+1; $hspp2=${$hsps{$besthit}}[$h2];$h2++) {
			if ($merge[$h2]) {next;}
			($s2,$e2)=split(/_/,$hspp2);
			my $frametest=abs($s1-$s2)/3;
if ($verbose) {    print "\tframetest: $frametest\n"; }
			next if ($frametest =~ /\D/);
#print "\t\t - trying to merge...\n";
			if ($s1<=$s2 && $e2<=$e1)	{ 
				$merge[$h2]="Y";
if ($verbose) {    print "\t\tmerging s1[s2e2]e1: ",$s1,"_",$e1," $hspp2 : ",$s1,"_",$e1,"\n"; }
			}
			elsif ($s2<=$s1 && $e1<=$e2)	{ 
if ($verbose) {    print "\t\tmerging s2[s1e1]e2: ",$s1,"_",$e1," $hspp2 : "; }
				$s1=$s2;
				$e1=$e2;
				$merge[$h2]="Y";
if ($verbose) {    print $s1,"_",$e1,"\n"; }
			}
			elsif ($s1<=$s2 && $s2<=$e1 && $e1<=$e2)	{ 
if ($verbose) {    print "\t\tmerging s1[s2e1]e2: ",$s1,"_",$e1," $hspp2 : "; }
				$e1=$e2;
if ($verbose) {     print $s1,"_",$e1,"\n"; }
				$merge[$h2]="Y";
			}
			elsif ($s2<=$s1 && $s1<=$e2 && $e2<=$e1)	{ 
if ($verbose) {    print "\t\tmerging s2[s1e2]e1 : ",$s1,"_",$e1," $hspp2 : "; }
				$s1=$s2;
if ($verbose) {    print $s1,"_",$e1,"\n"; }
				$merge[$h2]="Y";				
			}
			
			elsif ($linktrans==1 && $e1<=$s2) {
				my $link=$seqobj->subseq($e1,$s2);
				my $trans=1;
				for ($l=$e1;$l<$s2-2;$l+=3){
					my $codon=substr $link, $l, 3;
					if ($reverse) {
						if ($codon=~m/tca/i || $codon=~m/tta/i || $codon=~m/cta/i) { 
							$trans=0;
if ($verbose) {    	print "(rev) linkSTOP: $l : $codon\n";		}			
						}
					}
					else {
						if ($codon=~m/tga/i || $codon=~m/taa/i || $codon=~m/tag/i) { 
							$trans=0; 
if ($verbose) {    		print "linkSTOP: $l : $codon\n";		}			
						}
					}
				}
				if ($trans) {
if ($verbose) {    print "\t\tmerging over translatable in-frame gap: ",$s1,"_",$e1," $hspp2 : "; }
					$e1=$e2;
if ($verbose) {    print $s1,"_",$e1,"\n"; }
					$merge[$h2]="Y";				
				}
			}
		}
		push @rehsp, $s1."_".$e1; 
		push @recoord, $s1, $e1;
	}
#print "merging HSPs:\n@{$hsps{$hitname}}\n@rehsp\n";
#print "merging coord:\n@{$coord{$hitname}}\n@recoord\n";
	@{$hsps{$besthit}}=@rehsp;
	@{$coord{$besthit}}=@recoord;

#   print "\t@{$coord{$besthit}}\n";
#   print "\t@{$hsps{$besthit}}\n";

#----------------
# assembling CDS from HSPs and linkers

	next if (!$besthit);
   	my @coordinates=sort {$a <=> $b} @{$coord{$besthit}};
   	my @HSPs= sort {$a <=> $b} @{$hsps{$besthit}};
if ($verbose) {      	print "HSPs:\t\t@HSPs\n"; }
   	   	
   	my $ORF="";
   	my $h2=0;
   	my $skipquery;
if ($verbose) {    print "extracting:\t";	  	}
	for (my $h=0 ; $coordinates[$h] ; $h=$h+2) {
		
		my $linker=0;
		($start,$end)=split /_/, $HSPs[$h2];
#print "\t\tcoord: $coordinates[$h] - $coordinates[$h+1]\n";
#print "\t\thsp:   $start - $end\n";
		if ($end>$coordinates[$h+1]){
			$linker=3 * (1+ int(($end-$coordinates[$h+1])/3));
			$coordinates[$h+2]=$linker+$coordinates[$h+1] unless (!$coordinates[$h+2]);
			$coordinates[$h+1]=$end-$linker;
#print "\t\tadjst: $coordinates[$h] - $coordinates[$h+1]\n";
		}

		else {
			$linker=3 * (1+int(($coordinates[$h+2]-$coordinates[$h+1]-2)/3)) unless (!$coordinates[$h+2]);
		}
		$h2++;
		if ($linker < 0) { die "\n\nTROUBLE: negative linker: $linker\n";}
		if ($coordinates[$h+1]-$coordinates[$h]<=0){
if ($verbose) {  	print "\t\t\tapparent HSP confusion: zero lenght subsequence to be called; skipping this sequence\n"; }
			$skipquery=1;
		}
		last if ($skipquery);
if ($verbose) {    print "$coordinates[$h]_$coordinates[$h+1] "; }
if ($verbose) {    print "(linker ", $linker/3, " codons) " unless (!$linker);	 }	
		$ORF.=$seqobj->subseq($coordinates[$h],$coordinates[$h+1]);
		for (my $l=0;$l<$linker;$l++) { $ORF.="N";}
	}
	next if $skipquery;
if ($verbose) {    print "\n"; }

#------------
# checking for stops in the ORF, clipping

	if ($stopcheck){
		my @stops;
		
		for ($b=0;$b<length($ORF)-2;$b+=3) {
			my $codon=substr $ORF, $b, 3;
			if ($reverse) {
				if ($codon=~m/tca/i || $codon=~m/tta/i || $codon=~m/cta/i) { 
					push @stops, $b; 
if ($verbose) {    print "(rev) STOP: $b : $codon\n";	}				
				}
			}
			else {
				if ($codon=~m/tga/i || $codon=~m/taa/i || $codon=~m/tag/i) { 
					push @stops, $b; 
if ($verbose) {    print "STOP: $b : $codon\n";				}	
				}
			}
		}
						
		if ($stops[1]) {
if ($verbose) {    		print "\n\t\MANY STOPS! skipping query\n"; }
			next;
		}
		elsif ($stops[0] && $stops[0]<=length($ORF)/2) { 
			$ORF=substr $ORF, $stops[0]+3;
			($start,$end)=split /_/, $HSPs[0];
			my $newstart = $start+$stops[0]+3;
			if ($newstart>$end) {
if ($verbose) {    	print "\n\t\WHOLE 1st HSP is gone; skipping query\n"; }
				next;
			}
			else {
if ($verbose) {    print "first HSP:\t",$HSPs[0],"\n"; }
				$HSPs[0]=$newstart."_".$end;
if ($verbose) {    print "new first HSP:\t",$HSPs[0],"\n"; }
			}
		}
		elsif ($stops[0] && $stops[0]>length($ORF)/2) { 
			($start,$end)=split /_/, $HSPs[$#HSPs];
			my $newend = $end-(length($ORF)-$stops[0]);
			$ORF=substr $ORF, 0, $stops[0];
			if ($newend<$start) {
if ($verbose) {    	print "\n\t\WHOLE last HSP is gone; skipping query\n"; }
				next;
			}
			else {
if ($verbose) {    print "last HSP:\t",$HSPs[$#HSPs],"\n"; }
				$HSPs[$#HSPs]=$start."_".$newend;
if ($verbose) {    print "new last HSP:\t",$HSPs[$#HSPs],"\n"; }
			}
		}
	}

#------------
# extending ends

	if ($ends){
		($start,$endd)=split /_/, $HSPs[0];
		($startt,$end)=split /_/, $HSPs[$#HSPs];
#print "ends to extend: \n\t\tstart\t$start\n\t\tend\t$end\n";
		$cterm="";
		$nterm="";
		for ($e=1;;$e+=3){
			last if (($end+$e+2)>$seqobj->length );
			$codon=$seqobj->subseq($end+$e,$end+$e+2);
#print "C-cod: $codon\n";
			last if (length($codon)<3);
			if ($reverse){
				last if ($codon=~m/tca/i || $codon=~m/tta/i || $codon=~m/cta/i);
			}
			else { 
				last if ($codon=~m/tga/i || $codon=~m/taa/i || $codon=~m/tag/i);
			}
			$cterm.=$codon;
		}	
		for ($e=1;;$e+=3){
			last if (($start-$e-2)<1);
			$codon=$seqobj->subseq($start-$e-2,$start-$e);
#print "N-cod: $codon\n";
			last if (length($codon)<3);
			if ($reverse){
				last if ($codon=~m/tca/i || $codon=~m/tta/i || $codon=~m/cta/i);
			}
			else { 
				last if ($codon=~m/tga/i || $codon=~m/taa/i || $codon=~m/tag/i);
			}
			$codon.=$nterm;
			$nterm=$codon;
#print "$nterm\n";
		}	
if ($verbose) {    print "added 5':\t",(length($nterm)/3)," codons\n"; }
if ($verbose) {    print "added 3':\t",(length($cterm)/3)," codons\n"; }
		$ORFe=$nterm.$ORF.$cterm;
	}

#----------
# writing new sequence file, reverse-complementing if necessary, writing translations as a separate file

    my $CDS = Bio::Seq->new(-seq          => $ORF,
                     -description      => "extracted_CDS " . $result->query_description,
                     -display_id => $result->query_name,
                     -alphabet         => 'dna' );
	if ($reverse) { $CDS=$CDS->revcom;}
	my $PROT=$CDS->translate;	

    my $CDSe = Bio::Seq->new(-seq          => $ORFe,
                     -description      => "extracted_CDS " . $result->query_description,
                     -display_id => $result->query_name,
                     -alphabet         => 'dna' );
	if ($reverse) { $CDSe=$CDSe->revcom;}
	my $PROTe=$CDSe->translate;	
	
	$outdna->write_seq($CDS);
	$outpro->write_seq($PROT);
	$outdnaE->write_seq($CDSe);
	$outproE->write_seq($PROTe);

# writing hit length table
#	$besthit=~s/jgi\|Nemve1\|(\d+).+/$1/;
	print {HITS} $result->query_name,"\t",$qlength,"\t",length($ORF)/3,"\t",$besthit,"\t",$bh_length,"\n";

}
my $inxname=$ARGV[0].".index";
unlink $inxname;
