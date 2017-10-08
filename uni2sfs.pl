#!usr/bin/env perl

my $usage="

uni2sfs.pl : converts merged-uniqued 2bRAD data (output of mergeUniq.pl) to 1d and 2d SFS
             (since loci/SNPs remain undetermined the sfs should be analyzed as folded)

Arguments:

uni=[file]  : output of mergeUniq.pl
pops=[file] : two-column tab-delimited table of population designations for some or all
              samples in mergeUniq output. If omitted, all samples are considered to be 
              from same population (a single 1d SFS will be produced)
maxdp=100   : max depth in any one individual
minind=5    : min non-zero depth in any one individual    
mintot=5    : minimal total number of observation across all individuals    
sb=0.25     : strand bias cutoff: ratio of 2x reverse-complement reads to 
              total number of reads (NB! not Fisher test p-value!).

Output:

A series of files ready to import into dadi or moments:

pop1.sfs, pop2.sfs : 1d sfs 
pop1.pop2.sfs      : 2d sfs 

Mikhail V. Matz, matz\@utexas.edu

";

my $uni="";
my $pops="singlepop";
my $maxdp=100;
my $minind=5;
my $mintot=5;
my $sb=0.25;
if (" @ARGV "=~/uni=(\S+)/) { $uni=$1; } else { die $usage;}
if (" @ARGV "=~/pops=(\S+)/) { $pops=$1; } 
if (" @ARGV "=~/maxdp=(\d+)/) { $maxdp=$1;}
if (" @ARGV "=~/minind=(\d+)/) { $minind=$1;}
if (" @ARGV "=~/mintot=(\d+)/) { $mintot=$1;}
if (" @ARGV "=~/sb=(\S+)/) { $sb=$1;}

# reading ind2pops file
my %i2p={};
my %p2i={};
if ($pops ne "singlepop"){
	open POPS, $pops or die "\n\ncannot open pops file $pops\n\n";
	while (<POPS>) {
		chop;
		(my $in,my $po)=split("\t",$_);
		$i2p{$in}=$po;
		$p2i{$po}=$in;
	}
	my @inds=keys %i2p;
	my @poplist=keys %p2i;
	close POPS;
}
else {
	my @poplist=($pops);
}

open UNI, $uni or die "\n\ncannot open input file $uni\n\n";
my $first=1;
my %popindex={};
my %alls={};
my $nind=0;

my %counts ={};
my $badsb=0;
my $toohigh=0;
my $badindcover=0;
my $toofew=0;
my $monomorph=0;
my $goods=0;

while (<UNI>){
	chop;
	if ($first) {
		my @indnames=split("\t",$_);
		splice @indnames, 0, 3;
# recording line positions of individuals for each pop
		foreach my $p (@poplist) {
			for(my $i=0;$i<=$#indnames;$i++) {
				if ($i2p{$i}==$p || $p eq "singlepop") { 
					push @{$popindex{$p}},$i;
					push @alls,$i;
					$nind++;
				}
			}
		}
		$first++;
		next;
	}
	if ($nind==0) { die "\n\ndid not find any of the individuals listed in $pops\n";}
	else {
		my @data=split("\t",$_);
		my $total=$data[2];
		my $revcom=$data[3];
		splice @data, 0, 3;

# filtering
		if (2*$revcom/$total<$sb) {
			$badsb++;
			next;
		}
		if ($total/(1+$#data)>$maxdp){
			$toohigh++;
			next;
		}
		my $min=1e+4;
		my $sums=0;
		my $prod=1;
		foreach my $d (@data[@alls]) { 
			if ($d<$min && $d>0) { $min=$d;}
			$sums+=$d;
			$prod*=$d;
		}
		if ($min<$minind){
			$badindcover++;
			next;
		}
		if ($sums<$mintot){
			$toofew++;
			next;
		}
		if ($prod>0){
			$monomorph++;
			next;
		}
		$goods++;

# recording sfs entries
		for($p1=0;$p1<=$#poplist;$p1++) {
			my @pp1=@data[@{$popindex{$poplist[$p1]}}];
			for($p2=$p1;$p2<=$#poplist;$p2++) {
				my @pp2=@data[@{$popindex{$poplist[$p2]}}];
					my $p1c=0;
					foreach my $e (@pp1){
						if ($e>0) { $p1c++;}
					}
					if ($p1 != $p2){
						my $p2c=0;
						foreach my $e (@pp2){
							if ($e>0) { $p2c++;}
						}
						$sfs{$poplist[$p1]}{$poplist[$p2]}{$p1c}{$p2c}++;
					}
					else { $sfs{$poplist[$p1]}{$poplist[$p1]}{$p1c}{$p1c}++; }

					
			
		
			


