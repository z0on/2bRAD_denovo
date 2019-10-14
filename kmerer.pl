#!/usr/bin/perl

$usage= "

kmerer.pl  : 

subsets kmer fasta files output by jellyfish based on kmer counts.

prints to STDOUT

arguments:
1 : kmer.fa filename
2 : minimum number of occurences of kmer to keep it

Example:
kmerer.pl D5.kmers.fa 2

Mikhail Matz, matz\@utexas.edu			 
";

my $fq=shift or die $usage;
my $nn=shift or die $usage;;
my $name="";
my $name2="";
my $seq="";
my %seqs={};
open INP, $fq or die "cannot open file $fq\n";
while (<INP>) {
	if ($_=~/^>(\d+)\s/) {
		$name2=$1;
		if ($name>=$nn) {
			print ">$name\n$seq\n";
		}
		$seq="";
		$ll=0;
		$name=$name2;
	}
	else{
		chomp;
		$seq.=$_;
	}
}
$name2=$1; 
if ($name>=$nn) {
	print ">$name\n$seq\n";
}
