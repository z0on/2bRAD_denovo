#!/usr/bin/perl

my $usage= "

vcf2fakeGenome.pl: rewrites SNP coordinates in VCF according to the coordinates
in the concatenated genome file produced by concatFasta.pl

Arguments: 
vcf=[vcf file]
tab=[tab file made by concatFasta.pl]

" ;

if ("@ARGV"=~/tab=(\S+)/) { $tab=$1 } else { die $usage; } 
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1 } else { die $usage; }
open TAB, $tab or die "cannot find tab file $tab\n";
open VCF, $vcf or die "cannot find vcf file $vcf\n";

my %ch={};
my %co={};
my $l;
my $chrom;
my $coord;
my $coord2;

while (<TAB>) {
	($l, $chrom, $coord) = split("\t",$_);
	$ch{$l}=$chrom;
	$co{$l}=$coord;
}

my $fakename="f_".$vcf;
open OUT, ">$fakename" or die "cannot create $fakename\n";

my $l;
my $pos;

while (<VCF>) {
	if ($_=~/^#/) { print {OUT} $_ and next; }
	chop;
	@vline=split(/\t/,$_);
	if (!$ch{$vline[0]} | !$co{$vline[0]}) { warn "unrecognized locus: $vline[0]\n";} 
	$vline[1]=$vline[1]+$co{$vline[0]}-1;
	$vline[0]=$ch{$vline[0]};
	print {OUT} join("\t",@vline),"\n";
}
	
	

	
	
	
	