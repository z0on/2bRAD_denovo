#!/usr/bin/perl
my $usage="

Merges individual uniqued files produced by uniquerOne.pl

arg1: common extension of uniqued files 

minDP=[integer] minimum sequencing depth to consider a uniqie tag. Default 5.

Mikhail Matz, May 2014
matz\@utexas.edu

prints to STDOUT

";

# use Statistics::ChiSquare;

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

my $mincount=5;
if ("@ARGV"=~/minDP=(\d+)/) { $mincount=$1;}

$glob=shift or die $usage;
opendir THIS, ".";
@ins=grep /\.$glob$/,readdir THIS;

my %total={};
my %revcom={};
my %counts={};
my @samples=();
my $tag="";
my $tot="";
my $rc="";
my $ind="";
my $co="";
my $seen={};
my @tags=();

@ins=sort @ins;
@allins=@ins;

foreach $unifile (@ins) {
warn "\tfile $unifile\n";
$now=localtime;
warn "reading: $now\n";
	open INP, $unifile or die "cannot open file $unifile\n";
	my $frst=1;
	while (<INP>) {
		next if ($_=~/^seq/);
		chop;
		my $line=$_;
		($tag,$tot,$rc,$ind,$co)=split("\t",$line);
		$ind=~s/\..+//;
		if ($frst) { 
			push @samples,$ind;
			$frst=0;
		}
		if (!$seen{$tag}) {
			my $rc=rcom($tag);			
			$tag=$rc;
			$rc=$tot-$rc;
		}
		if (!$seen{$tag}) {
			push @tags, $tag;
			$seen{$tag}=1;
		}
		$counts{$tag,$ind}=$tot;
		$revcom{$tag}+=$rc;
		$total{$tag}+=$tot;
	}
}
close INP;

$now=localtime;
warn "processing: $now\n";

sub bycount {
	$total{$b} <=> $total{$a}
}
@tags=sort bycount @tags;

open FAS, ">mergedUniqTags.fasta" or die "cannot create file mergedUniqTags.fasta";
print "tag\tseq\tcount\trevcom\t",join("\t",@samples),"\n";
my $globalcount=0;
my $indexx=0;
foreach $tag (@tags) {
	next if ($tag=~/HASH/ || $tag!~/[ATGC]+/);
	last if ($total{$tag}<$mincount);
	$globalcount+=$total{$tag};
	$indexx++;
	my $fahead="tag".$indexx." count=".$total{$tag};
	print {FAS} ">$fahead\n$tag\n";
	print "tag".$indexx."\t$tag\t$total{$tag}\t$revcom{$tag}";
	foreach $ind (@samples) {
		if ($counts{$tag,$ind}) { print "\t$counts{$tag,$ind}"; }
		else { print "\t0";}
	}
	print "\n";
}
close FAS;
warn "
considering reads seen at least $mincount times:
$globalcount reads processed
";
$now=localtime;
warn "$now\n";