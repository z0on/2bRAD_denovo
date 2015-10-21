#!/usr/bin/perl

my $usage="
removeBayescanOutliers.pl 

deletes (or extracts) rows from the vcf file corrsponding to outliers in 
bayescan _fst.txt output file

arguments:

bayescan=[filename]   : BayeScan *_fst.txt file name
vcf=[filename]        : vcf file name
FDR=[float]           : false discovery rate cutoff (q-value in BayeScan file), default 0.1
mode=[delete|extract] : whether to delete (default) or extract the outlier rows

";

my $vcf;
my $bs;
my $fdr=0.1;
my $mode="delete";

if ("@ARGV"=~/bayescan=(\S+)/) { $bs=$1;}
else {die $usage;}
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else {die $usage;}
if ("@ARGV"=~/FDR=(\S+)/) { $fdr=$1;}
if ("@ARGV"=~/mode=extract/) { $mode="extract";}

open BS, $bs or die "cannot open bayescan file $bs\n";
open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @splits=();
my $qval;
my @outs=();
while(<BS>) {
	chomp;
	@splits=split(/\s/,$_);
	$qval=$splits[4];
# warn "q: $qval\n";
	if ($qval<$fdr) { 
		push @outs, $splits[0];
		warn "outlier : ",$_,"\n";
	}
}
close BS;
if ($#outs==0) { exit "no outliers found\n";}

my $count=0;
while (<VCF>) {
	if ($_=~/^#/) { 
		print $_; 
		next;
	}
	$count++;
	if (" @outs "=~/ $count /) {
		if ($mode eq "delete") {
			warn "$count : ",$_;
		}
		else {
			print $_;
		}
	}
	else {
		if ($mode eq "delete") {
			print $_;
		}
	}
}		
