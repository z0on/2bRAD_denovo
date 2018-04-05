#!/usr/bin/perl

my $usage= "

ktable2loci.pl :

assembles locus-annotated table from mergeKmers.pl and cd-hit results
reverse-complements to align all kmers belonging to the same locus
filters by read depth, minInd and minIndCover

arg1: mergeKmers.pl result
arg2: cd-hit .clstr file
arg3: cd-hit .fasta file

          minDP = 2 : minimum total depth (count for the kmer in the whole dataset). 
         minInd = 2 : minumum number of individuals the kmer must be present in. 
    minIndCover = 1 : minimum kmer representation in an individual, if >0.
  mismatchCheck = 0 : check clustered kmers for excessive mismatches?
mismatchLimit = 0.9 : warn if clustered kmers show more than this fraction of mismatches.

prints to STDOUT

Example: 
ktable2loci.pl kmers.ktab cdh_kmers.fasta.clstr cdh_kmers.fasta >mydata.ltab

";

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

sub mism16 {
	my @ref=split("",$_[0]);
	my @test=split("",$_[1]);
	my $mism=0;
	my $i=0;
	for ($i=0;$i<16;$i++){
		if ($test[$i] ne $ref[$i]) { $mism++;}
	}
	return $mism;
}

sub mism {
	my @ref=split("",$_[0]);
	my @test=split("",$_[1]);
	my $mism=0;
	my $i=0;
	for ($i=0;$i<=$#ref;$i++){
		if ($test[$i] ne $ref[$i]) { $mism++;}
	}
	$mism=$mism/($i+1);
	return $mism;
}


my $ufile=shift or die $usage;
my $cfile=shift or die $usage;
my $fafile=shift or die $usage;
my $count=2;
my $minind=2;
my $minindcov=1;
my $mismatchCheck=0;
my $mismatchLimit=0.9;
if ("@ARGV"=~/minInd=(\d+)/) { $minind=$1;}
if ("@ARGV"=~/minIndCover=(\d+)/) { $minindcov=$1;}
if ("@ARGV"=~/minDP=(\d+)/) { $count=$1;}
if ("@ARGV"=~/mismatchLimit=(\d+)/) { $mismatchLimit=$1;}
if ("@ARGV"=~/mismatchCheck=(\d)/) { $mismatchCheck=$1;}

open FA, $fafile or die "cannot open clusters file $cfile\n";
my %tag2seq={};
my $tag="";
my $seq="";
while (<FA>) {
	chomp;
	if ($_=~/^>(\S+)/){
		if ($seq) {	$tag2seq{$tag}=$seq;}
		$tag=$1;
		$seq="";
	}
	else { $seq.=$_;}
}
$tag2seq{$tag}=$seq;
close FA;

open CL, $cfile or die "cannot open clusters file $cfile\n";
my %tag2cl={};
my %ref={};
my $cl;
my $first=0;
while (<CL>){
	chop;
	if ($_=~/^>Cluster (\d+)/) {
#warn "Cluster $1 ...\n";
		$cl=$1;
		$first=1;
	}
	elsif ($_=~/\s>(k\d+)/) {
		$tag2cl{$1}=$cl;
		if ($first){
			$ref{$cl}=$tag2seq{$1};
			$first=0;
		}
	}
}
undef %tag2seq;
close CL;

open UNI, $ufile or die "cannot open kmers table $ufile\n";

my $tag;
my $tagseq;
my $tagcount;
my @rest;
my $locus;
my $dpcount=0;
my $minindcount=0;
my $minindcovcount=0;
my $total=0;
my %vars={};

while(<UNI>){
	if ($_=~/kmer/) {
		print "locus\t",$_;
		next;
	}
	chop;
	$total++;
	($tag,$tagseq,$tagcount,@rest)=split("\t",$_);
	if (exists($tag2cl{$tag}) && $ref{$tag2cl{$tag}}) { 
		$locus=$tag2cl{$tag};
		my $refseq=$ref{$tag2cl{$tag}};
		my $mins=1e+6;
		my $nz=0;
		for (my $i=0;$i<=$#rest;$i++) {
			if ($rest[$i]>0) { 
				$nz++;
				if ($rest[$i]<$mins) { $mins=$rest[$i];}
			}
		}
		if ($tagcount<$count) {	next; }
		$dpcount++;
		if($nz<$minind){ next;}
		$minindcount++;
		if($mins<$minindcov) { next;}
		$minindcovcount++;
		my $mi=mism16($refseq,$tagseq);
		if ($mi>6) { 
#warn "revcom?\nref $refseq\ntst $tagseq\nmism $mi\n";
			$tagseq=rcom($tagseq);
			$revcom=$tagcount-$revcom;
		}
		if ($mismatchCheck!=0) {
			$mi=mism($refseq,$tagseq);
			if ($mi>$mismatchLimit) { 
				warn "\tpossibly incompatible kmers in cluster $tag2cl{$tag}:\n\tref\t$refseq\n\ttest\t$tagseq\n\tmism $mi\n\n";
				next;
			}				
		}
		push @{$vars{$locus}}, "locus$locus\t$tag\t$tagseq\t$tagcount\t",join("\t",@rest),"\n";
	}
	else { die "\nstrage $tag: seq:$tag2seq{$tag} cluster:$tag2cl{$tag}  ref:$ref{$tag2cl{$tag}}\n";}
}	
close UNI;
	
foreach my $locus (sort {$a <=> $b} keys %vars){
	next if ($locus=~/HASH/);
	foreach my $var (@{$vars{$locus}}){
		print $var;
	}
}

print STDERR "
dataset $ufile:
$total\ttotal tags
$dpcount\tpass depth cutoff $count
$minindcount\tpass minInd cutoff $minind
$minindcovcount\tpass minIndCover cutoff $minindcov

";
	
	
	