#!/usr/bin/perl
$usage="

Trims 2b-RAD fastq files to leave only those with  adaptor on the far end and restriction site

arg1: fastq filename
arg2: pattern defining the site, such as \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\' for BcgI
arg3: (optional) adaptor sequence, default AGATCGGAAGA
arg4: (optional) number of bases to trim off the ends of the reads (corresponding to ligated overhangs), default 0

";

open INP, $ARGV[0] or die $usage;
my $site=$ARGV[1] or die $usage;
my $adap="AGATCGGAAGA";
if ($ARGV[2]){$adap=$ARGV[2];}
my $clip=0;
if ($ARGV[3]){$clip=$ARGV[3];}
my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";
my $ll=3;
while (<INP>) {
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		if ($seq=~/^($site)$adap/) {
#			if (!$trim) {$trim=length($1)-$clip*2;}
			my $rd=substr($1,$clip,length($1)-$clip*2);
			print "$name\n$rd\n+\n",substr($qua,$clip,length($1)-$clip*2),"\n";
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
$name2=$1; 
if ($seq=~/^($site)$adap/) {
	if (!$trim) {$trim=length($1);}
	print "$name\n$1\n+\n",substr($qua,0,$trim),"\n";
}
