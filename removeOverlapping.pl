#!/usr/bin/env perl

my $in="";
$in=shift @ARGV;
open INP, $in or die "cannot open input file (arg1) $in\n";

my $line1=<INP>;
chomp $line1;
print $line1,"\n";
my $line2="";
while (<INP>) {
	chomp;
	my $line2=$_;
	(my $chra,my $a1, my $a2)=split('[\:\-]',$line1);
	(my $chrb,my $b1, my $b2)=split('[\:\-]',$line2);
	if ($chra ne $chrb ) { 
		print $line2,"\n";
		$line1=$line2;
		next; 
	}
	if ($b1 > $a2) { 
		print $line2,"\n";
		$line1=$line2;
		next; 
	}
	warn "------\noverlap:\n",$line1,"\n",$line2,"\n";
}		
		
	
	