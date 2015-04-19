#!/usr/bin/perl

my $usage="

thinner.pl : 

Leaves only one SNP per given interval, the one with the highest minor allele frequency.

Arguments:

vcf=[file name]    : vcf file name
interval=[integer] : interval length, default 40 (for 2bRAD tags)

Output:  thinned VCF, printed to STDOUT

Example: thinner.pl vcf=denovo.recal.vcf > thinDenov.vcf

NOTE: do not thin your variants if you want to calculate Tajima's D!
Also, to select variants for dadi analysis, use --thin option in vcftools instead.

";

my $vcf;
if (" @ARGV "=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage; }
my $inter=40;
if (" @ARGV "=~/interval=(\d+)/) { $inter=$1; } 

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @dats;
my $info;
my $chrom;
my $pos;
my $pos0=1;
my $toout;
my $af;
my $maxaf=0;

while (<VCF>) {
	if ($_=~/^#/) { print and next;}
	chomp;
	@dats=split(/\t/,$_);
	$chrom=$dats[0];
	$pos=$dats[1];
	$info=$dats[7];
	if ($info=~/AF=(\S*?);/) { $af=$1; } else { warn "no AF:\n@dats\n" and next;}
	if ($af>0.5) { $af=1-$af; }
	if ($pos-$pos0 <=$inter and ($chrom eq $chrom0)){
		if ($af >$maxaf) {
			$maxaf=$af;
			$toout=join("\t",@dats);
		}
		$pos0=$pos;
	}
	else {
		print "$toout\n" unless (!$toout);
		$pos0=$pos;
		$chrom0=$chrom;
		$maxaf=$af;
		$toout=join("\t",@dats);
	}
}
print "$toout\n" unless (!$toout);
