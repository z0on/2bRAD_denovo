#!/usr/bin/perl
my $usage="

uniquerOne.pl :

Makes uniqued 2bRAD read file for a single fastq file.
This is analogous to making 'stacks' in STACKS.
The script records unique tag sequences, the total
number of their appearances, and how many of those
were in reverse-complement orientation.
(STACKS would consider reverse-complements as separate loci)

prints to STDOUT

arg1: fastq filename

Example: 
uniquerOne.pl G4_TCAG.trim >G4_TCAG.trim.uni

";

$infile=shift or die $usage;
my @ins=($infile);

my %count={};
my %inds={};
my %qual={};
my %seen={};
my %seenall={};
my %revs;

sub rcom {
	my $rcs=scalar reverse ("$_[0]");
	$rcs=lc $rcs;
	$rcs=~s/a/T/g;
	$rcs=~s/t/A/g;
	$rcs=~s/g/C/g;
	$rcs=~s/c/G/g;
#	$rcs=~s/![atgcATGC]/N/g;
	return $rcs;
}

my $now;

foreach $gmap (@ins) {
warn "\tfile $gmap\n";
$now=localtime;
warn "$now\n";
	open INP, $gmap or die "cannot open file $gmap\n";
#	$gmap=~s/\.$glob//;
	my $ii=0;
	my $qu="";
	my $seq="";
	my $rc="";
	my $rev=0;
	my $ii=1;
	my $ccc=0;
	while (<INP>) {
		chop;
		my $line=$_;
		$ccc++;
		if ($ii==1) {
			if (!$qu) {
				$ii=2;
#warn "no quality for seq|$seq|: $qu , line $ccc\n";
				next;
			}
#warn "\t\tquality for seq|$seq|: $qu , line $ccc\n";
			if (!$seenall{$seq}) {
				$rc=rcom($seq);
				if ($seenall{$rc}){ 
					$seq=$rc;
					$rev=1;
					$qu=scalar reverse ("$qu");
				}
				else { $seenall{$seq}=1;}
			}
			$count{$seq}{$gmap}++;
			push @{$qual{$seq}}, $qu;
			push @{$inds{$seq}}, $gmap unless $seen{$seq}{$gmap};
			push @{$revs{$seq}},$rev;
#print "\n$seq\n$qu\t\t$count{$seq}{$gmap}\n";	
			$seen{$seq}{$gmap}=1;
			$qu="";
			$seq="";
			$rev=0;
			$ii=2;
		}
		elsif ($ii==2) {
			$seq=$line;
			$ii=3;
		}
		elsif ($ii==3){ 
			$ii=4;
		}
		else {
			$qu=$line;
			$ii=1;
		}
	}
}
if (!$seenall{$seq}) {
	$rc=rcom($seq);
	if ($seenall{$rc}){ 
		$seq=$rc;
		$rev=1;
		$qu=scalar reverse ("$qu");
	}
	else { $seenall{$seq}=1;}
}
$count{$seq}{$gmap}++;
push @{$qual{$seq}}, $qu;
push @{$inds{$seq}}, $gmap unless $seen{$seq}{$gmap};
push @{$revs{$seq}},$rev;
#print "\n$seq\n$qu\t\t$count{$seq}{$gmap}\n";	
$seen{$seq}{$gmap}=1;

$now=localtime;
warn "processing: $now\n";

sub bycount {
	$#{$qual{$b}} <=> $#{$qual{$a}}
}

my @seqs=keys %inds;
@seqs= sort bycount @seqs;

my $globalcount=0;
print "seq\ttotal\trevcom\tsamples\tcounts\n";
foreach $seq (@seqs) {
	next if ($seq=~/HASH/);
#	my @avg=();
#	foreach $q ( @{$qual{$seq}}) {
#		my @qc=split("",$q);
#		for($i=0;$i<=$#qc;$i++){
#			$avg[$i]+=(ord("$qc[$i]")-33);
#		}
#print "$q\n";
#	}
#	for($i=0;$i<=$#avg;$i++){ 
#		$avg[$i]=$avg[$i]/($#{$qual{$seq}}+1);
#		$avg[$i]=sprintf("%.0f",$avg[$i]/($#{$qual{$seq}}+1));
#		$avg[$i]=chr(($avg[$i]<=93? $avg[$i] : 93) + 33);
#	}
#print join("",@avg),"\n\n";
	my @counts;
	foreach $ind (@{$inds{$seq}}){
		push @counts, $count{$seq}{$ind};
	}
	my $rsum=0;
	foreach $r (@{$revs{$seq}}){ $rsum+=$r;}
#print "$rsum\t",($#{$qual{$}}+1-$rsum),"\t";
#	my $chi=chisquare($rsum,($#{$qual{$seq}}+1-$rsum));
#	if ($chi=~/<(\d+)/) { $chi=$1;}
#print $chi,"\n";
	$globalcount+=($#{$qual{$seq}}+1);
	print $seq,"\t",($#{$qual{$seq}}+1),"\t$rsum\t",
	join(",",@{$inds{$seq}}),"\t",join(",",@counts),"\n";
}
warn "\n$globalcount reads processed\n";
$now=localtime;
warn "$now\n";
