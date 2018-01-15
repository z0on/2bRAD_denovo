#!/usr/bin/env perl
my $usage="

realsfs2dadi.pl :

Coverts output of 'realSFS dadi' command to dadi-snp format understood by dadi and moments.

argument 1: input file generated with
realSFS dadi pop1.saf.idx pop2.saf.idx -sfs pop1.sfs -sfs pop2.sfs -ref reference.fasta -anc ancestral.fasta >dadiout

subsequent arguments: numbers of individuals in each population (must match the number and order of populations 
analyzed by realSFS)

prints to STDOUT

Example:
realsfs2dadi.pl dadiout 20 20 >mypops_dadi.data

Mikhail Matz, matz\@utexas.edu
";

my $infile=shift @ARGV or die $usage;
if ($#ARGV==0) { die $usage; }
my @popsizes=@ARGV;

open IN, $infile or die "cannot open input file $infile\n\n $usage";
print "REF\tOUT\tAllele1";
for (my $i=0;$i<=$#popsizes;$i++) { print "\tpop",$i; }
print "\tAllele2";
for (my $i=0;$i<=$#popsizes;$i++) { print "\tpop",$i; }
print "\tGene\tPosition\n";

while (<IN>) {
	chomp;
	(my $ref, my $anc, my $chr, my $pos, my @popdata)=split("\t",$_);
	my @ra=split("",$ref);
	my @counts=();
	for ($i=0; $i<$#popdata;$i+=2){	$counts[$i/2]=$popdata[$i]; }
	print "$ref\t$anc\t$ra[1]";	
	for (my $i=0;$i<=$#popsizes;$i++) { print "\t",2*$popsizes[$i]-$counts[$i]; }
# printing fake derived allele:
	if ($ra[1]=~/[AGT]/){ print "\tC";} else { print "\tT";}
	for (my $i=0;$i<=$#popsizes;$i++) { print "\t",$counts[$i]; }
	print "\t$chr\t$pos\n";
}
