#!/usr/bin/perl

my $usage="

counts number of SNPs per tag position in a vcf file made by haplocall_denovo

arguments:

vcf=[file name]  input vcf file

";

my $vcf;
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage;}

open VCF, $vcf or die "cannot open input file $vcf\n";

my %snpcount={};
while (<VCF>){
	if ($_=~/TP=(\d+)/) {$snpcount{$1}++;}
}
for (my $i=0;$i<36;$i++){
	print "$i\t$snpcount{$i}\n";
}