#!/usr/bin/perl

my $usage="

filter_stats.pl :

Calculates quantiles for three INFO fields in a VCF file:
- total coverage depth (DP)
- strand bias (SB)
- allele bias (AB)

For AB and SB, low values are bad and high values are good;
for DP, both very low and very high values are bad, the best values
should be around the median.

Arguments:

vcf=[file name]  vcf file to analyze

Output: 
Tables of quantiles (printed to STDOUT).

Example:
filter_stats.pl vcf=alltags.ul_Variants_count10_ab10_sb10_clip0.vcf

NOTE: only use this to select your filtering criteria if you don't have 
genotyping replicates.
 
";

my $vcf;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @dp=();
my @ab=();
my @sb=();
my @dline=();
my @iline=();
my $geno;

while (<VCF>) {
	next if ($_=~/^#/);
	@dline=split("\t",$_);
	$geno=()=$_=~/\d\/\d/g;
	next if ($geno==0);
	@iline=split(";",$dline[7]);
	foreach my $i (@iline){
		if ($i=~/AB=(\d+)/) { push @ab, $1;}
		elsif ($i=~/SB=(\d+)/) { push @sb, $1;}
		elsif ($i=~/DP=(\d+)/) { push @dp, sprintf("%.0f",$1/$geno);}
	}
}
close VCF;

@ab=sort {$a <=> $b} @ab;
@sb=sort {$a <=> $b} @sb;
@dp=sort {$a <=> $b} @dp;

my $total=$#dp+1;

my %abq;
my %sbq;
my %dpq;
my %dpq2;

$abq{1}=$ab[sprintf("%.0f",$total*0.01)];
$sbq{1}=$sb[sprintf("%.0f",$total*0.01)];
$abq{5}=$ab[sprintf("%.0f",$total*0.05)];
$sbq{5}=$sb[sprintf("%.0f",$total*0.05)];
$abq{10}=$ab[sprintf("%.0f",$total*0.1)];
$sbq{10}=$sb[sprintf("%.0f",$total*0.1)];
$abq{15}=$ab[sprintf("%.0f",$total*0.15)];
$sbq{15}=$sb[sprintf("%.0f",$total*0.15)];
$abq{20}=$ab[sprintf("%.0f",$total*0.2)];
$sbq{20}=$sb[sprintf("%.0f",$total*0.2)];
$abq{30}=$ab[sprintf("%.0f",$total*0.3)];
$sbq{30}=$sb[sprintf("%.0f",$total*0.3)];
$abq{40}=$ab[sprintf("%.0f",$total*0.4)];
$sbq{40}=$sb[sprintf("%.0f",$total*0.4)];
$abq{50}=$ab[sprintf("%.0f",$total*0.5)];
$sbq{50}=$sb[sprintf("%.0f",$total*0.5)];
$abq{60}=$ab[sprintf("%.0f",$total*0.6)];
$sbq{60}=$sb[sprintf("%.0f",$total*0.6)];
$abq{70}=$ab[sprintf("%.0f",$total*0.7)];
$sbq{70}=$sb[sprintf("%.0f",$total*0.7)];
$abq{80}=$ab[sprintf("%.0f",$total*0.8)];
$sbq{80}=$sb[sprintf("%.0f",$total*0.8)];
$abq{90}=$ab[sprintf("%.0f",$total*0.9)];
$sbq{90}=$sb[sprintf("%.0f",$total*0.9)];


$dpq{1}=$dp[sprintf("%.0f",$total*0.005)];
$dpq2{1}=$dp[sprintf("%.0f",$total*0.995)];
$dpq{5}=$dp[sprintf("%.0f",$total*0.025)];
$dpq2{5}=$dp[sprintf("%.0f",$total*0.975)];
$dpq{10}=$dp[sprintf("%.0f",$total*0.05)];
$dpq2{10}=$dp[sprintf("%.0f",$total*0.95)];
$dpq{15}=$dp[sprintf("%.0f",$total*0.075)];
$dpq2{15}=$dp[sprintf("%.0f",$total*0.925)];
$dpq{20}=$dp[sprintf("%.0f",$total*0.1)];
$dpq2{20}=$dp[sprintf("%.0f",$total*0.9)];
$dpq{30}=$dp[sprintf("%.0f",$total*0.15)];
$dpq2{30}=$dp[sprintf("%.0f",$total*0.85)];
$dpq{40}=$dp[sprintf("%.0f",$total*0.2)];
$dpq2{40}=$dp[sprintf("%.0f",$total*0.8)];
$dpq{50}=$dp[sprintf("%.0f",$total*0.25)];
$dpq2{50}=$dp[sprintf("%.0f",$total*0.75)];
$dpq{60}=$dp[sprintf("%.0f",$total*0.3)];
$dpq2{60}=$dp[sprintf("%.0f",$total*0.7)];
$dpq{70}=$dp[sprintf("%.0f",$total*0.35)];
$dpq2{70}=$dp[sprintf("%.0f",$total*0.65)];
$dpq{80}=$dp[sprintf("%.0f",$total*0.4)];
$dpq2{80}=$dp[sprintf("%.0f",$total*0.6)];
$dpq{90}=$dp[sprintf("%.0f",$total*0.45)];
$dpq2{90}=$dp[sprintf("%.0f",$total*0.55)];

print "\nAB quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %abq){
	print 100-$q,"%:\t>$abq{$q}\n";
}

print  "\nSB quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %sbq){
	print 100-$q,"%:\t>$sbq{$q}\n";
}

print  "\nmeanDP quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %dpq){
	print  100-$q,"%:\t$dpq{$q} - $dpq2{$q}\n";
}


