#!/usr/bin/perl

my $usage="

dephase.pl: leaves a single SNP out of a phase set (for STRUCTURE)

arg1    input vcf
ps=3    position of the Phase Set identifier in the genotype string

prints to STDOUT

";

my $pspos=3;
if ("@ARGV"=~/ps=\d+/) { $pspos=$1;}

open VCF, $ARGV[0] or die "specify VCF file to dephase\n";

my @dline=();
my %seen={};
while (<VCF>) {
	if ($_=~/^#/) { 
		print $_;
		next;
	}
	@dline=split("\t",$_);
	my @gt=split(":",$dline[9]);
	my $ps=$gt[$pspos-1];
	if ($ps ne "." && $seen{$ps}) { next; }
	else {
		print $_;
		$seen{$ps}=1;
	}
}