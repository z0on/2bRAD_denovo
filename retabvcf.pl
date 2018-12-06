#!/usr/bin/perl

my $usage ="

retabvcf.pl :

rewrites coordinates of variants in a VCF file according to a table
of locations of original contigs in pseudo-chromosomes
(created by concatFasta.pl)

Arguments:

vcf=[file name] : vcf file to re-coordinate
tab=[file name] : tab-delimited table of contig locations in the form
                  contig-chromosome-start-end

Output: 
VCF with reformatted coordinates (printed to STDOUT)

Example: 
retabvcf.pl vcf=gatk_after_vqsr.vcf tab=where/genome/is/mygenome_cc.tab > retab.vcf

Mikhail Matz, matz\@utexas.edu, October 2014

";

my $vcf="";
my $tab;
if (" @ARGV "=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage; }
if (" @ARGV "=~/tab=(\S+)/) { $tab=$1; } else { die $usage; }

my %ctg={};
my %start={};
my $ctig="";
my $chr="";
my $beg="";
my $end="";
open TAB, $tab or die "cannot open table $tab\n";
while(<TAB>){
	chop;
	($ctig,$chr,$beg,$end)=split(/\t/,$_);
	push @{$starts{$chr}},$beg;
	$ctg{$chr}{$beg}=$ctig;
}
close TAB;

open VCF, $vcf or die "cannot open vcf file $vcf\n";
my @vline=();
while (<VCF>) {
	if ($_=~/^#/) { 
		print $_ ;
		next;
	}
	chop;
	@vline=split(/\t/,$_);
	$chr=$vline[0];
	$beg=$vline[1];
	my $i=0;
	++$i until ${$starts{$chr}}[$i] > $beg or $i > $#{$starts{$chr}};
	my $cbeg=${$starts{$chr}}[$i-1];
	$vline[0]=$ctg{$chr}{$cbeg};
	$vline[1]=$beg-$cbeg+1;
#warn "$chr $beg : $vline[0] $vline[1]\n";
	print join("\t",@vline),"\n";
}