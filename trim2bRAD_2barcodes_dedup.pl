#!/usr/bin/perl
$usage="

trim2bRAD_2barcodes_dedup.pl : 

This script does three things:

- Filters 2bRAD fastq reads to leave only those with a 100% matching restriction site,
degenerate 5'-leader, secondary 3'-barcode, and adaptor on the far 3'end, 
trims away everything except the IIb-fragment itself;

- Deduplicates: removes all but one read sharing the same 64-fold degenerate
leader, the first 34 bases of the insert sequence, and secondary barcode
(this results in 128-fold dynamic range: 64-fold degeneracy x 2 strand orientations);

- Splits reads into separate files according to secondary barcode.

Writes trimmed fastq files named according to the secondary barcodes detected. 

Arguments:
input=[fastq filename]
site=[perl-style pattern for restrictase site] 
		Default is BcgI: \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\'
bc2=[perls-style barcode pattern] in-line barcode that immediately follows the RAD fragment, 
		default \'[ATGC]{4}\'
adaptor=[far-end adaptor sequence] default AGATC
clip=[integer] number of bases to clip off the ends of the reads, default 0
minBCcount=[integer] minimum count per in-line barcode to output a separate file. 
				     Default 100000.
sampleID=[integer] the position of name-deriving string in the file name
					if separated by underscores, such as: 
					for input file Sample_RNA_2DVH_L002_R1.cat.fastq
					specifying arg2 as \'3\' would create output 
					file with a name \'2DVH.trim'

Example:
trim2bRAD_2barcodes_dedup.pl input=Myproject_L8_G3.fastq sampleID=3 

";

my $input="";
if ("@ARGV"=~/input=(\S+)/){ $input=$1; }
else { die $usage; }
open INP, $1 or die "cannot open input $1\n";
my $site='.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}';
if ("@ARGV"=~/site=(\S+)/){ $site=$1;}
my $bcod="[ATGC]{4}";
if ("@ARGV"=~/bc2=(\S+)/){ $bcod=$1;}
my $degen="[ATGC][ATGC][AG][AT]CC";
my $adap="AGATC";
if ("@ARGV"=~/adaptor=(\S+)/){ $adap=$1;}
my $clip=0;
if ("@ARGV"=~/clip=(\d+)/){ $clip=$1;}
my $sampleid=100;
if("@ARGV"=~/sampleID=(\d+)/){ $sampleid=$1;}
my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";
my $minBCcount=100000;
if ("@ARGV"=~/minBCcount=(\d+)/){$minBCcount=$1;}

my $ll=3;
my %data={};
my %dups={};
my $counter=0;
my %dupnum={};
my $chk;
my $bcc;
my $read;
while (<INP>) {
	chomp;
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		$counter++;
#print "$seq:";
		if ($seq=~/^($degen)($site)($bcod)$adap/ and $seq!~/N/) {
			$read=$2;
			$bcc=$3;
			my $chk1=substr($read,0,34);
			$chk=$1.$chk1.$bcc;
#print "$1:$2:$3|check:$chk\n";
			if ($dups{$chk}) {
				$dupnum{$bcc}++;
#print "\t\tdup found: $chk\n";
			}
			else {
				$dups{$chk}=1;
				my $rd=substr($read,$clip,length($read)-$clip*2);
				$qua=substr($qua,$clip,length($read)-$clip*2);
				$dline="$name bcd=$bcc\n$rd\n+\n$qua\n";
				push @{$data{$bcc}}, $dline ;
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
                        $bcc=$3; 
                        my $chk1=substr($read,0,34);
                        $chk=$1.$chk1.$bcc;
#print "$1:$2:$3|check:$chk\n";
                        if ($dups{$chk}) {
                                $dupnum{$bcc}++;
#print "\t\tdup found: $chk\n";
                        }
                        else {
                                $dups{$chk}=1;
                                my $rd=substr($read,$clip,length($read)-$clip*2);
                                $qua=substr($qua,$clip,length($read)-$clip*2);
                                $dline="$name bcd=$bcc\n$rd\n+\n$qua\n";
                                push @{$data{$bcc}}, $dline ;
                        }
#print "$dline\n";
                }
                else {
#print "\n\nNOT FINDING $site:$bcod:$adap\n\n";
}


my $outname;
if ($sampleid<100) {
	my @parts=split('_',$input);
	$outname=$parts[$sampleid-1];
}
else { 
	$input=~s/\..+//;
	$outname=$input;
}

my $totgood=0;
my $totdups=0;
foreach my $bc2 (sort keys %data) {
	next if ($bc2=~/HASH/);
	next if ($#{$data{$bc2}}<$minBCcount);
	print $outname,"_","$bc2: goods: ",$#{$data{$bc2}}+1," ; dups: $dupnum{$bc2}\n";
	$totgood=$totgood+$#{$data{$bc2}}+1;
	$totdups=$totdups+$dupnum{$bc2};
	my $outname1=$outname."_".$bc2.".tr0";
	open OUT, ">$outname1" or die "cannot create $outname1\n";
	foreach $d (@{$data{$bc2}}) { print {OUT} $d; }
	close OUT;
}

print "\n$outname: total goods : $totgood ; total dups : $totdups\n";