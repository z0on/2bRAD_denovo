#!/usr/bin/perl

my $usage="

thinner.pl (v.2) : 

Leaves only one SNP per given interval, either the one with the highest minor allele 
frequency, or chosen at randomom.

Arguments:

   vcf=[file name] : vcf file name
interval=[integer] : interval length, default 40 (for 2bRAD tags)
    randomom=[0|1] : whether to choose SNPs randomomly. With 0 (default), the SNP with 
                     the highest minor allele frequency will be chosen.

Output:  thinned VCF, printed to STDOUT

Example: thinner.pl vcf=denovo.recal.vcf > thinDenov.vcf

NOTES: - Do not thin variants if you want to calculate pi or Tajima's D.
       - For dadi, use with option random=1 

Mikhail Matz, matz@utexas.edu, 09/16/2015

";

my $vcf;
if (" @ARGV "=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage; }
my $inter=40;
if (" @ARGV "=~/interval=(\d+)/) { $inter=$1; } 
my $random=0;
if (" @ARGV "=~/random=1/) { $random=1; } 

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @dats;
my $info;
my $chrom;
my $pos;
my $pos0=1;
my $toout;
my $af;
my $maxaf=0;
my @snps=();
my $chrom0="";

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
		push @snps,join("\t",@dats);
		if ($af >$maxaf) {
			$maxaf=$af;
			$toout=join("\t",@dats);
		}
	}
	else {
		if ($random) { print $snps[rand @snps],"\n" unless (!$snps[0]); }
		else { print "$toout\n" unless (!$toout); }
		$pos0=$pos;
		$chrom0=$chrom;
		$maxaf=$af;
		$toout=join("\t",@dats);
		@snps=();
		push @snps,join("\t",@dats);
	}
}
if ($random) { print $snps[rand @snps] unless (!$snps[0]); }
else { print "$toout\n" unless (!$toout); }
