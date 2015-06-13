#!/usr/bin/perl
my $usage="

mergeUniq.pl :

Merges individual uniqued files produced by uniquerOne.pl into one table.

arg1: common extension of uniqued files 
minDP=[integer] minimum sequencing depth to consider a uniqie tag. Default 5.

prints to STDOUT

Example:
mergeUniq.pl uni >mydataMerged.uniq

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

my $tag;
my $tot1;
my $rc1;
my $ind1;
my $tot2;
my $rc2;
my $chi1;
my $chi2;
my $counts1;
my $counts2;
my $ind2;
my %seen={};
my $now;
my %dstring={};
my @rest;

my @ins=sort @ins;
my @allins=();

my $frst=1;

$reff=shift(@ins);
open INP, $reff or die "cannot open first file $reff\n"; 
warn "\tfile $reff\n";
$now=localtime;
warn "reading: $now\n";
while (<INP>) {
	next if ($_=~/^seq/);
	chop;
	($tag,@rest)=split("\t",$_);
	@{$dstring{$tag}}=@rest;
	if ($frst) {
		push @allins,$rest[2];
		$frst=0;
	}
}

foreach $gmap (@ins) {
	warn "\tfile $gmap\n";
	$frst=1;
	$now=localtime;
	warn "reading: $now\n";
	open INP, $gmap or die "cannot open file $gmap\n";
	while (<INP>) {
		next if ($_=~/^seq/);
		chop;
		my $line=$_;
		($tag,$tot2,$rc2,$ind2,$counts2)=split("\t",$line);
		if ($frst) {
			push @allins,$ind2;
			$frst=0;
		}
		if (!$dstring{$tag}) {
			my $rc=rcom($tag);			
			if ($dstring{$rc}){ 
				$tag=$rc;
				$rc2=$tot2-$rc2;
			}
		}
		if (!$dstring{$tag}) { @{$dstring{$tag}}=($tot2,$rc2,$ind2,$counts2);}
		else {
			$tot1=$tot2+${$dstring{$tag}}[0];
			$rc1=$rc2+${$dstring{$tag}}[1];
#warn "$tag|@{$dstring{$tag}}\n";
#			$chi1=chisquare($rc1,($tot1-$rc1));
#warn "$chi1\n\n";
#			if ($chi1=~/<(\d+)/) { $chi1=$1;}
#			else { die "chisquare error: $chi\n";}
			$ind1=${$dstring{$tag}}[2].",".$ind2;
			$counts1=${$dstring{$tag}}[3].",".$counts2;
			@{$dstring{$tag}}=($tot1,$rc1,$ind1,$counts1);
		}
	}
}
#warn "ALL allins: @allins\n";
$now=localtime;
warn "processing: $now\n";

sub bycount {
	${$dstring{$b}}[0] <=> ${$dstring{$a}}[0]
}

my @seqs=keys %dstring;
@seqs= sort bycount @seqs;

foreach $s (@allins){
	$s=~s/\..+//;
}

#@allins

open FAS, ">mergedUniqTags.fasta" or die "cannot create file mergedUniqTags.fasta";
print "tag\tseq\tcount\trevcom\t",join("\t",@allins),"\n";
my $globalcount=0;
my $indexx=0;
foreach $seq (@seqs) {
	next if ($seq=~/HASH/ || $seq!~/[ATGC]+/);
	last if (${$dstring{$seq}}[0]<$mincount);
	$globalcount+=${$dstring{$seq}}[0];
	$indexx++;
	my $fahead="tag".$indexx." count=".${$dstring{$seq}}[0];
	print FAS ">$fahead\n$seq\n";
	my @cts=split(",",pop @{$dstring{$seq}});
	my @inds=split(",", pop @{$dstring{$seq}});
#warn "$seq inds: @inds\n";
#warn "$seq counts: @cts\n";
	my @ctall=();
	my $skip=0;
	for ($i=0;$ind=$allins[$i];$i++) {
		$inds[$i-$skip]=~s/\..+//;
#warn "\ti $i, skip $skip. is ind $inds[$i-$skip] the same as $ind?\n";
		if ($inds[$i-$skip] eq $ind) { 
#warn "\tyes\n\n";
			push @ctall, $cts[$i-$skip];
		}
		else {
			push @ctall, 0;
			$skip++;
#warn "\tno\n";
		}
	}
#warn "\t\tcounts: @ctall\n\n";
	push @{$dstring{$seq}}, @ctall;
	print "tag$indexx\t$seq\t",join("\t",@{$dstring{$seq}}),"\n";
	undef $dstring{$seq};
}
close FAS;
warn "\nconsidering reads seen at least $mincount times:\n$globalcount reads processed\n";
$now=localtime;
warn "$now\n";
