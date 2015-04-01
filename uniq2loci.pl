#!/usr/bin/perl

my $usage= "

uniq2loci.pl :

assembles locus-annotated table from uniquer and cd-hit results
filters by read depth and strand bias (to avoid memory issues in genotype calling)
flips orientation of cluster members to match the most abundant tag in a cluster

saves cluster members seen just once and their references (singl.tab file) to use for 
error model development (requires previous running of mergeUniq.p with minDP=1)

arg1: uniquer file
arg2: cd-hit .clstr file
arg3: cd-hit .fasta file

minDP=[integer] minimum read depth (total counts for the tag in the whole dataset). 
                Default: 5
minSB=[integer] strand bias cutoff ( min(direct:reverse,reverse:direct)x100 ). 
                Default: 5 

prints to STDOUT

Example: 
uniq2loci.pl mydataMerged.uniq cdh_alltags.fas.clstr cdh_alltags.fas >cdh_alltags.ul

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

sub mism12 {
	my @ref=split("",$_[0]);
	my @test=split("",$_[1]);
	my $mism=0;
	for (my $i=0;$i<12;$i++){
		if ($test[$i] ne $ref[$i]) { $mism++;}
	}
	return $mism;
}

my $ufile=shift or die $usage;
my $cfile=shift or die $usage;
my $fafile=shift or die $usage;
my $skew=5;
my $count=5;
if ("@ARGV"=~/minSB=(\d+)/) { $skew=$1;}
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

open UNI, $ufile or die "cannot open uniquer file $ufile\n";

my $tag;
my $tagseq;
my $tagcount;
my $revcom;
my @rest;
my $locus;
my $dpcount=0;
my $sbcount=0;
my $total=0;
my %vars={};

open SI, ">singl.tab" or die "cannot create singl.tab\n";

while(<UNI>){
	if ($_=~/revcom/) {
		print "locus\t",$_;
		next;
	}
	chop;
	$total++;
	($tag,$tagseq,$tagcount,$revcom,@rest)=split("\t",$_);
	if (exists($tag2cl{$tag}) && $ref{$tag2cl{$tag}}) { 
		$locus=$tag2cl{$tag};
		my $refseq=$ref{$tag2cl{$tag}};
		my $sb=$revcom/($tagcount-$revcom);
		if ($sb>1){ $sb=1/$sb;}
		$sb=sprintf("%.0f",$sb*100);
		if ($tagcount<$count) {
			if ($tagcount==1 && $refseq ne $tagseq) { print {SI} "$refseq\t$tagseq\n";}
			next;
		}
		$dpcount++;
		next if ($sb<$skew);
		$sbcount++;
		my $mi=mism12($refseq,$tagseq);
		if ($mi>4) { 
#warn "revcom?\nref $refseq\ntst $tagseq\nmism $mi\n";
			$tagseq=rcom($tagseq);
			$revcom=$tagcount-$revcom;
		}
		$mi=mism12($refseq,$tagseq);
		if ($mi>4) { 
#			warn "\tincompatible tags:\n\tref\t$refseq\n\ttest\t$tagseq\n\tmism $mi\n";
			next;
		}				
		push @{$vars{$locus}}, "locus$locus\t$tag\t$tagseq\t$tagcount\t$revcom\t",join("\t",@rest),"\n";
	}
	else { die "\nstrage $tag: seq:$tag2seq{$tag} cluster:$tag2cl{$tag}  ref:$ref{$tag2cl{$tag}}\n";}
}	
close SI;
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
$sbcount\tpass strand bias cutoff $skew
";
	
	
	