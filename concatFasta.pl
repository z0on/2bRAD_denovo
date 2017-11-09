#!/usr/bin/perl


my $usage="

concatFasta.pl : 

concatenates fasta records in a draft genome reference into \"pseudo-chromosomes\"
(either a specified number of them or as many as needed to 
split the genome into chunks just over 100G bases in length) 

Do this to your genotyping reference to improve the speed and memory usage 
of the UnifiedGenotyper (GATK)

Arguments:

fasta=[file name] : fasta file name
num=[integer]     : number of pseudo-chromosomes, default 10

Output:

[fasta filename base]_cc.fasta : new fasta file
[fasta filename base]_cc.tab   : table of original sequence IDs, their coordinates
                                 in pseudo-chromosomes, and full original header lines

Example:
concatFasta.pl fasta=mygenome.fasta

";

my $fas="";
if (" @ARGV "=~/fasta=(\S+)/) { $fas=$1 } else { die $usage; }
my $nc = 10;
if (" @ARGV "=~/num=(\d+)/) { $nc=$1; }
#my $bl = 100;
#if (" @ARGV "=~/blocks=(\d+)/) { $bl=$1; }
my $maxlen=100000000;

my $l=`grep ">" $fas | wc -l`;
chop $l;
my $chunk=sprintf("%.0f",$l/$nc);

my $namebase=$fas;
$namebase=~s/\..+//;
my $outname=$namebase."_cc.fasta";
my $tabname=$namebase."_cc.tab";
open OUTF, ">$outname" or die "cannot create $outname\n";
open OUTT, ">$tabname" or die "cannot create $outname\n";

if ($chunk*$nc<$l) { $chunk++;}
print "\nconcatenating $l records into $nc pseudo-chromosomes 
$chunk records per chromosome\n"; 
my $k=1;
my $chrom="chr".$k;

open INP, $fas or die "cannot open input $fas\n";
my $counter=0;

my $seq="";
my $seq1="";
my $title="";
my $defline="";

while (<INP>) {
	chomp;
	if ($_=~/^>/ ) { 
		if ($seq1) {
			$seq1=~s/\s//;
			print {OUTT} "\t",length($seq)+1,"\t",length($seq)+length($seq1),"\t",$defline,"\n";
			$seq.=$seq1;
			$counter++;
			if (length($seq)>$maxlen) { goto DUMP; }
#warn "$chrom\t",length($seq),"\n";
		}
		if ($counter>$chunk) {
			DUMP:
			$k++;
			$seq=~s/(.{100})/$1\n/g; # wraps lines in blocks of 100
			print {OUTF} ">$chrom\n$seq\n";
			$seq="";
			$counter=1;
			$chrom="chr".$k;
		}
		$seq1="";
		my @tit;
		@tit=split(/\s+/,$_);
		my $title=shift @tit;
		$title=~s/>//;
		my $defline=join(" ",@tit);
		print {OUTT} "$title\t$chrom";
	}
	else { 	$seq1.=$_; 	}
}
$seq1=~s/\s//;
print {OUTT} "\t",length($seq)+1,"\t",length($seq)+length($seq1),"\t",$defline,"\n";
$seq=~s/(.{100})/$1\n/g;
print {OUTF} ">$chrom\n$seq\n";
