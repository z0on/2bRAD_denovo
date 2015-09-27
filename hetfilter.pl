#!/usr/bin/perl

my $usage="
hetfilter.pl :

Filters out loci with too many heterozygotes, to guard against lumped paralogous loci.
Also filters out loci with too many missing genotypes.

Arguments:

vcf=[file name]  input vcf file
maxmiss=[float]  allowed fraction of missing genotypes, default 0.5.
 maxhet=[float]  maximum fraction of heterozygotes among non-missing genotypes.
                 Default 0.75
                 
Prints to STDOUT

Mikhail Matz, matz\@utexas.edu

";

my $vcf;
my $maxmiss=0.5;
my $maxhet=0.75;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/maxmiss=(\S+)/) { $maxmiss=$1;}
if ("@ARGV"=~/maxhet=(\S+)/) { $maxhet=$1;}

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my $badHeteros=0;
my $badMissing=0;
my $totals=0;
my $goods=0;

while (<VCF>) {
	if ($_=~/^#/) { print and next;}
	chop;
	$totals++;
	my @lin=split("\t",$_);
	my @start=splice(@lin,0,9);
	my $Missing = () = "@lin" =~ /\.[\/\|]\./gi;
	my $Heteros = () = "@lin" =~ /0[\/\|]1/gi;
	my $Refhomos = () = "@lin" =~ /0[\/\|]0/gi;
	my $Nsam= scalar @lin;
	if ($Missing/$Nsam > $maxmiss) { 
warn "\tHet:$Heteros;Miss:$Missing;Tot:$Nsam - badMiss\n";
		$badMissing++;
		next;
	}
	if ($Heteros/($Nsam-$Missing+0.01) > $maxhet) { 	
warn "Het:$Heteros;Miss:$Missing;Tot:$Nsam - badHet\n";
		$badHeteros++;
		next; 	
	}
	$goods++;
	print join("\t",@start)."\t".join("\t",@lin)."\n";
}
warn "
$totals total loci
$badMissing dropped because fraction of missing genotypes exceeded $maxmiss
$badHeteros dropped because fraction of heterozygotes exceeded $maxhet
$goods written

";