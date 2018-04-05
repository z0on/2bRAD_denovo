#!/usr/bin/perl

my $usage= "

ktable2loci.pl :

assembles locus-annotated table from mergeKmers.pl and cd-hit results
filters by read depth, minInd and minIndCover

arg1: mergeKmers.pl result
arg2: cd-hit .clstr file
arg3: cd-hit .fasta file

minDP=[integer] minimum read depth (total counts for the tag in the whole dataset). 
                Default: 5
minInd=[integer] number of individuals the tag must be present in. 
                Default: 2
minIndCover=[integer] minimum representation in an individual (if >0). Default 1.

prints to STDOUT

Example: 
ktable2loci.pl mydata.ktab cdh_kmers.fasta.clstr cdh_kmers.fasta >mydata.ltab

";


my $ufile=shift or die $usage;
my $cfile=shift or die $usage;
my $fafile=shift or die $usage;
my $count=5;
my $minind=2;
my $minindcov=1;
if ("@ARGV"=~/minInd=(\d+)/) { $minind=$1;}
if ("@ARGV"=~/minIndCover=(\d+)/) { $minindcov=$1;}
if ("@ARGV"=~/minDP=(\d+)/) { $count=$1;}

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
	elsif ($_=~/\s>(tag\d+)/) {
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
	if ($_=~/revcom/) {
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
		for (my $i=0;$rest[$i];$i++) {
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
		my $mi=mism12($refseq,$tagseq);
		if ($mi>4) { 
#			warn "\tincompatible tags:\n\tref\t$refseq\n\ttest\t$tagseq\n\tmism $mi\n";
			next;
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
$total\ttotal tags
$dpcount\tpass depth cutoff $count
$minindcount\tpass minInd cutoff $minind
$minindcovcount\tpass minIndCover cutoff $minindcov
";
	
	
	