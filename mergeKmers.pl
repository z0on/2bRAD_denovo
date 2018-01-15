#!/usr/bin/env perl
my $usage="

mergeKmers.pl :

Generates counts table from multiple kmer fasta files produced by Jellyfish.

To go from fastq to kmer counts (for one file):
fastq_to_fasta -i mydata1.fastq -o mydata1.fasta
jellyfish count -m [kmer length] -s 100M -t [number of threads] -C mydata1.fasta -o mydata1.jf
jellyfish dump mydata1.jf > mydata1_kmers.fa

arg1: list of kmer fasta files
minDP=[integer] minimum total sequencing depth to consider a uniqie kmer. Default 5.
minInd=[integer] minimum number of individuals the kmer must be seen in. Default 2.
maxInd=[integer] maximum number of individuals the kmer must be seen in. Default: 
                 (total number of individuals - minInd). Set this to twice per-pop sample 
                 size for kmer sharing analysis.

prints to STDOUT

Example:
mergeKmers.pl kmerfiles minDP=5 minInd=5 maxInd=50 > kmers_5to50.txt

";

#sub rcom {
#	my $rcs=scalar reverse ("$_[0]");
#	$rcs=lc $rcs;
#	$rcs=~s/a/T/g;
#	$rcs=~s/t/A/g;
#	$rcs=~s/g/C/g;
#	$rcs=~s/c/G/g;
##	$rcs=~s/![atgcATGC]/N/g;
#	return $rcs;
#}

my $mincount=5;
my $minind=2;
if ("@ARGV"=~/minDP=(\d+)/) { $mincount=$1;}
if ("@ARGV"=~/minInd=(\d+)/) { $minind=$1;}

my $flist=shift or die $usage;
open FL, $flist or die "cannot open file list $flist\n";
my @ins;
while (<FL>){
	chomp;
	push @ins, $_;
}

my $maxind=$#ins+1-$minind;
if ("@ARGV"=~/maxInd=(\d+)/) { $maxind=$1;}

my $tag;
my $tot;
my $ind;
my $counts;
my $now;
my %dstring={};

my @allins=sort @ins;

$reff=shift(@ins);
open INP, $reff or die "cannot open file $reff\n"; 
$now=localtime;
warn "reading input files...  $now\n";
warn "$reff\n";
while (<INP>) {
	chomp;
	if ($_=~/^>(\d+)/) {
		$counts=$1;
	}
	else { 
		$tag=$_;
		@{$dstring{$tag}}=($counts,$reff,$counts);
	}
}

foreach $gmap (@ins) {
	warn "$gmap\n";
	$frst=1;
	$now=localtime;
	open INP, $gmap or die "cannot open file $gmap\n";
	while (<INP>) {
		chomp;
		if ($_=~/^>(\d+)/) {
			$counts=$1;
		}
		else { 
			$tag=$_;
#			if (!$dstring{$tag}) {
#				my $rc=rcom($tag);			
#				if ($dstring{$rc}){ 
#					$tag=$rc;
#				}
#			}
			if (!$dstring{$tag}) { @{$dstring{$tag}}=($counts,$gmap,$counts);}
			else {
				$tot=$counts+${$dstring{$tag}}[0];
				$ind=${$dstring{$tag}}[1].",".$gmap;
				$counts=${$dstring{$tag}}[2].",".$counts;
				@{$dstring{$tag}}=($tot,$ind,$counts);
			}
		}
	}
}

$now=localtime;
warn "processing: $now\n";

sub bycount {
	${$dstring{$b}}[0] <=> ${$dstring{$a}}[0]
}

my @seqs=keys %dstring;
@seqs= sort bycount @seqs;


my $globalcount=0;
my $indexx=0;
print "id\tkmer\tcount\t",join("\t",@allins),"\n";
foreach $seq (@seqs) {
	next if ($seq=~/HASH/ || $seq!~/[ATGC]+/);
	last if (${$dstring{$seq}}[0]<$mincount);
	$globalcount+=${$dstring{$seq}}[0];
	$indexx++;
	my @cts=split(",",pop @{$dstring{$seq}});
	my @inds=split(",", pop @{$dstring{$seq}});
	if ($#inds+1<$minind) {	next; }
	elsif ($#inds+1>$maxind) { next;}	
	my @ctall=();
	my $skip=0;
	for ($i=0;$ind=$allins[$i];$i++) {
		if ($inds[$i-$skip] eq $ind) { 
			push @ctall, $cts[$i-$skip];
		}
		else {
			push @ctall, 0;
			$skip++;
		}
	}
	push @{$dstring{$seq}}, @ctall;
	print "k$indexx\t$seq\t",join("\t",@{$dstring{$seq}}),"\n";
	undef $dstring{$seq};
}
warn "\nconsidering kmers seen at least $mincount times:\n$globalcount kmer occurences\n";
$now=localtime;
warn "done   $now\n\n";
