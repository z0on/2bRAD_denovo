#!/usr/bin/perl

my $usage = "

alleleCounts.pl: 

Uses vcf file to count alleles for each locus across populations.
Can filter counted genotypes by genotype quality.
Only counts bialelic loci.

Arguments:

 vcf=[file name] : input vcf v.4 file
pops=[file name] : tab-delimited table: sample name <tab> population affiliation
 minGQ=[integer] : minimum phred-scaled genotype quality to count alleles. Default 0.

prints to STDOUT:
locus number
chromosome
position
ref allele
alt allele
population
ref allele count
alt allele count
number of missing alleles

Mikhail V. Matz, matz/@utexas.edu

";

my $mingq=0;
my $vcf;
my $pops;
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage;}
if ("@ARGV"=~/pops=(\S+)/) { $pops=$1; } else { die $usage;}
if ("@ARGV"=~/minGQ=(\d+)/) { $mingq=$1; } 

open VCF, "$vcf" or die "cannot open input vcf file $vcf\n";
open POPS, "$pops" or die "cannot open input populations file $pops\n";

my %inds2pops={};
my %pops2inds={};

while (<POPS>) {
	chop;
	(my $ind,my $pop)=split('\t',$_);
	$inds2pops{$ind}=$pop;
	push @{$pops2inds{$pop}}, $ind;
}
my @Pops = keys %pops2inds ;
my $GQpos=-1;
my %num2pops={};
my %pops2num={};
my $totals=0;
print "SNP\tCHROM\tPOS\tREF\tALT\tPOP\tREFC\tALTC\tMISSC\n";

while (<VCF>) {
	if ($_=~/^#CHROM/) { 
		chop;
		my @lin=split("\t",$_);
		my @start=splice(@lin,0,9);
		for(my $i=0;$lin[$i];$i++) {
			if ($inds2pops{$lin[$i]}) {
				$num2pops{$i}=$inds2pops{$lin[$i]};
				push @{$pops2num{$inds2pops{$lin[$i]}}}, $i;
			}
			else { warn "\t no pop affiliation found for $lin[$i]\n";}
		}
		next;
	}

	if ($_=~/^#/) { next;}
	
	chop;
	my @lin=split("\t",$_);
	my @start=splice(@lin,0,9);
	if ($start[4]=~/\,/) { next;} # multi-allelic SNP
	if ($GQpos==-1 && $mingq>0) {
#warn "format: $start[8]\n";
		if ($start[8]=~/GQ/) {
			my @info=split('\:',$start[8]);
			my $i;
			for ($i=0;$info[$i];$i++){ 
				if ($info[$i] eq "GQ") { last;}
			}
			$GQpos=$i;
#warn "GQpos:$GQpos\n";
		}
		else { 
			warn "\tno GQ field found\n"; 
			$mingq=0;
		}
	}
	my %popRef={};
	my %popAlt={};
	my %popMiss={};
	$totals++;

	for($i=0;$lin[$i];$i++){
		my $pop = $num2pops{$i};
#warn "$i:pop $pop; geno $lin[$i]\n";
		my @geno=split("\:",$lin[$i]);
		if (($mingq>0 && $geno[$GQpos]<$mingq) || $geno[0] eq "./." ) {
				$popMiss{$pop}+=2;
		}
		else {
			if ($geno[0] eq "0/0" || $geno[0] eq "0|0" ) { 
				$popRef{$pop}+=2;
			}
			elsif ($geno[0] eq "1/1" || $geno[0] eq "1|1" ) { 
				$popAlt{$pop}+=2;
			}
			elsif ($geno[0] eq "0/1" || $geno[0] eq "0|1" ) { 
				$popRef{$pop}++;
				$popAlt{$pop}++;
			}
			elsif ($geno[0] eq "1/0" || $geno[0] eq "1|0" ) { 
				$popRef{$pop}++;
				$popAlt{$pop}++;
			}
			else { warn "weird genotype: $geno[0]\n"; }
		}
	}

	foreach my $p (@Pops) {
		next if ($p=~/HASH/);
		if (!$popRef{$p}){ $popRef{$p}=0;}
		if (!$popAlt{$p}){ $popAlt{$p}=0;}
		if (!$popMiss{$p}){ $popMiss{$p}=0;}
		print "$totals\t$start[0]\t$start[1]\t$start[3]\t$start[4]\t$p\t$popRef{$p}\t$popAlt{$p}\t$popMiss{$p}\n";
	}
}
warn "total $totals loci\n";