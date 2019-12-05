#!/usr/bin/perl
$usage="

trim2bRAD_dedup.pl : 

This is a reduced version of main trimmer script that does not split data by secondary
(in-read) barcode. 

- Filters 2bRAD fastq reads to leave only those with a 100% matching restriction site,
degenerate 5'-leader, secondary 3'-barcode, and adaptor on the far 3'end, 
trims away everything except the IIb-fragment itself;

- Deduplicates: removes all but one read sharing the same 64-fold degenerate
leader, the first 34 bases of the insert sequence, and secondary barcode
(this results in 128-fold dynamic range: 64-fold degeneracy x 2 strand orientations);

Writes trimmed reads to STDOUT. 

Arguments:
input=[fastq filename]
site=[perl-style pattern for restrictase site] 
		Default is BcgI: \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\'
adaptor=[far-end adaptor sequence] default AGATC?
deduplicate=1  whether to remove PCR duplicates; set to 0 to keep them.
clip=[integer] number of bases to clip off the ends of the reads, default 0

Example:
trim2bRAD_dedup.pl input=Myproject_L8_G3.fastq > Myproject_L8_G3.tr0

";

my $input="";
if ("@ARGV"=~/input=(\S+)/){ $input=$1; }
else { die $usage; }
open INP, $1 or die "cannot open input $1\n";
my $site='.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}';
if ("@ARGV"=~/site=(\S+)/){ $site=$1;}
my $degen="[ATGC][ATGC][AG][AT]CC";
my $adap="AGATC?";
if ("@ARGV"=~/adaptor=(\S+)/){ $adap=$1;}
my $clip=0;
if ("@ARGV"=~/clip=(\d+)/){ $clip=$1;}
my $sampleid=100;
if("@ARGV"=~/sampleID=(\d+)/){ $sampleid=$1;}
my $deduplicate=1;
if ("@ARGV"=~/deduplicate=0/){ $deduplicate=0;}
my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";

my $ll=3;
my %data={};
my %dups={};
my $counter=0;
my %dupnum={};
my $chk;
my $read;
while (<INP>) {
	chomp;
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		$counter++;
#print "$seq:";
		if ($seq=~/^($degen)($site)$adap/ and $seq!~/N/) {
			$read=$2;
			if ($deduplicate==0) { $read=$1.$2;} 
			my $chk1=substr($read,0,34);
			$chk=$1.$chk1;
#print "$1:$2:$3|check:$chk\n";
			if ($dups{$chk} and $deduplicate) {
				$dupnum++;
#print "\t\tdup found: $chk\n";
			}
			else {
				$dups{$chk}=1;
				my $rd=substr($read,$clip,length($read)-$clip*2);
				$qua=substr($qua,$clip,length($read)-$clip*2);
				$dline="$name\n$rd\n+\n$qua\n";
				push @{$data}, $dline ;
			}
#print "$dline\n";
		}
		else {
#print "\n\nNOT FINDING $site:$bcod:$adap\n\n";
		}
#print "-----------\n$counter\n$_\n";
		$seq="";
		$ll=0;
		$qua="";
		@sites=();
		$name=$name2;
	}
	elsif ($ll==0){
		chomp;
		$seq=$_;
		$ll=1;
	}
	elsif ($ll==2) { 
		chomp;
		$qua=substr($_,3,length($_)); 
		$ll=3;
	}
	else { $ll=2;}
}

if ($seq=~/^($degen)($site)($bcod)$adap/ and $seq!~/N/) {
                        $read=$2;
                        my $chk1=substr($read,0,34);
                        $chk=$1.$chk1;
#print "$1:$2:$3|check:$chk\n";
                        if ($dups{$chk}) {
                                $dupnum++;
#print "\t\tdup found: $chk\n";
                        }
                        else {
                                $dups{$chk}=1;
                                my $rd=substr($read,$clip,length($read)-$clip*2);
                                $qua=substr($qua,$clip,length($read)-$clip*2);
                                $dline="$name\n$rd\n+\n$qua\n";
                                push @{$data}, $dline ;
                        }
#print "$dline\n";
                }
                else {
#print "\n\nNOT FINDING $site:$bcod:$adap\n\n";
}


my $totgood=0;
my $totdups=0;
$totgood=$totgood+$#{$data}+1;
$totdups=$totdups+$dupnum;
foreach $d (@{$data}) { print $d; }

warn "\n$input: goods:$totgood dups:$totdups\n";