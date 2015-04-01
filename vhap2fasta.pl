#!/usr/bin/perl

my $usage= "

vhap2fasta.pl :

extracts reference 2bRAD tags from a vcf \"Vhap\" file 
(produced by haplocall_denovo.pl) in fasta format

Arguments:
arg1: vhap.vcf file name
altalleles=[0|1] whether to output alternative alleles. Default 0.

prints to STDOUT

";

my $vcf=shift or die $usage;
my $altalleles=0;
if ("@ARGV"=~/altalleles=1/) { $altalleles=1;}

open INP, $vcf or die "cannot open input file $vcf\n";

my $locus;
my $pos;
my $major;
my $dot;
my $alts;
my $info;
my @rest;

while (<INP>) {
	chop;
	next if ($_=~/^#/);
	@rest=split("\t",$_);
	$locus=$rest[0];
	$pos=$rest[1];
	$major=$rest[3];
	@alts=split(",",$rest[4]);
	$info=$rest[7];
	if ($info=~/AF=(\S+)$|AF=(\S+);/) { @afs=split(":",$1);}
	print ">",$locus,"_",$pos,"_maj AF=",$afs[0],"\n$major\n";
	if ($altalleles==1){
		next if ($alts[0] eq ".");
		for(my $i=0;$i<=$#alts;$i++){
			print ">",$locus,"_",$pos,"_alt",$i," AF=",$afs[$i+1],"\n$alts[$i]\n";
		}
	}
}
