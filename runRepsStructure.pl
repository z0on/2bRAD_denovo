#!/usr/bin/perl

my $usage = "
runRepsStructure.pl :
makes a command file for running 100 replicate runs of fastStructure

Argument line must be just like for fastStructure, for example
(assuming the input file is myStrucFile.str and it is STRUCTURE format)

runRepsStructure.pl --input myStrucFile --output myStr -K 5 --format str --prior logistic

";

my $fname="";
if ("@ARGV"=~/output\s+(\S+)/) { $fname=$1;}
else { die $usage;}

my $argline0="python \$FAST_STRUCTURE_DIR/sucture.py @ARGV";
for ($i=1;$i<=100;$i++) {
	my $subout=$fname.$i;
	my $argline=$argline0;
	$argline=~s/\s$fname\s/ $subout /;
	print $argline,"\n";
}
