#!/usr/bin/perl
my $usage="

removes reads containing N calls from a fastq file
arg1: fastq filename
prints to STDOUT

Mikhail Matz, matz\@utexas.edu, September 2014

";

open INP, $ARGV[0] or die $usage;
my $name="";
my $name2="";
my $seq="";
my $qua="";
my $ll=3;
while (<INP>) {
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		if ($seq and $seq!~/N/) {
			print "$name\n$seq\n+\n$qua\n";
		}
		$seq="";
		$ll=0;
		$qua="";
		@sites=();
		$name=$name2;
	}
	elsif ($ll==0){
		chomp;
		$seq=$_;
		$ll=1;
	}
	elsif ($ll==2) { 
		chomp;
		$qua=$_;
		$ll=3; 
	}
	else { $ll=2;}
}
if ($seq and $seq!~/N/) {
	print "$name\n$seq\n+\n$qua\n";
}
