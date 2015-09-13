#!/usr/bin/perl

my $usage="

replicatesMatch.pl : (version 0.6, September 13 2015)

Selects polymorphic variants that have identical genotypes among replicates.  

Genotypes of two replicates should be identical or missing; the maximal fraction of 
genotype pairs with missing genotypes is controlled by missing=  argument.

Arguments:

vcf=[file name]  input vcf file
replicates=[file name] - a two column tab-delimited table listing 
                 pairs of samples that are replicates (at least 3 pairs)

matching=[float] required fraction of matching genotypes 
                 (missing counts as match if there are other non-missing matching pairs), 
                 default 1 (all must match)

missing=[float]  allowed fraction of missing genotypes, default 0.25

altPairs=2       minimal number of matching genotypes involving non-reference alleles. 

hetPairs=0       minimal number of matching heterozygotes. 

max.het=0.5      maximum fraction of heterozygotes among non-missing genotypes
                 (guards against lumped paralogous loci).

polyonly=[1|0]   extract only polymorphic sites. Default 0.

Example: 
replicatesMatch.pl vcf=round2.vcf replicates=clonepairs.tab > vqsr.vcf

";

my $vcf;
my $reps;
my $missing=0.25;
my $fmatch=1;
#my $allalts=0;
my $maxhet=0.5;
my $hetPairs=0;
my $altPairs=2;
my $polyonly=0;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/replicates=(\S+)/) { $reps=$1;}
else { die $usage; }
if ("@ARGV"=~/missing=(\S+)/) { $missing=$1;}
if ("@ARGV"=~/matching=(\S+)/) { $fmatch=$1;}
#if ("@ARGV"=~/allAlts=1/ ) { $allalts=1;}
if ("@ARGV"=~/max.het=(\S+)/) { $maxhet=$1;}
if ("@ARGV"=~/hetPairs=(\d+)/) { $hetPairs=$1;}
if ("@ARGV"=~/altPairs=(\d+)/) { $altPairs=$1;}
if ("@ARGV"=~/polyonly=1/) { $polyonly=1;}

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @samples;
my @pairs=();
my %indr={};
my @npairs=();
my $r1;
my $r2;
my $nreps=0;
my $pass=0;
my $hetpass=0;
my $altpass=0;
my $total=0;
my $poly=0;
my $numalt=0;
my $numout=0;
my $numhets=0;

while (<VCF>) {
	if ($_=~/^#/) {
		if ($_=~/contig/) { next;}
		elsif ($_=~/^#CHROM/){
			print $_;
			chop;
			@samples=split("\t",$_);
			my @lead=splice (@samples,0,9);
			open REP, $reps or die "cannot open replicates file $reps\n";
			while (<REP>){
				next if ($_!~/\S+/);
				$nreps++;
				chomp;
				($r1,$r2)=split(/\s+/,$_);
				my $collect=0;
				for(my $i=0;my $s=$samples[$i];$i++){
					if ($s eq $r1) { 
						$indr{$r1}=$i;
						$collect++;
					}
					elsif ($s eq $r2) { 
						$indr{$r2}=$i;
						$collect++;
					}
				}
				push @pairs,"$indr{$r1}:$indr{$r2}";
				push @npairs,"$r1:$r2";
				if ($collect<2) { die "cannot locate samples $r1 and/or $r2\n";}
			}
			close REP;
			$minfalt=1/$nreps unless $nreps==0;
#warn "$nreps replicate pairs\n";
#warn join("\t",@npairs),"\n",join("\t",@pairs),"\n";
		}
		else { print $_;}
		next;
	}
	chop;
	$total++;
	my @lin=split("\t",$_);
	my @start=splice(@lin,0,9);
	
	my $Missing = () = "@lin" =~ /\.[\/\|]\./gi;
#	my $Heteros = () = "@lin" =~ /0[\/\|]1/gi;
#	my $Althomos = () = "@lin" =~ /1[\/\|]1/gi;
#	my $Refhomos = () = "@lin" =~ /0[\/\|]0/gi;
	my $Nsam= scalar @lin;	
	if ($Heteros/($Nsam-$Missing) > $maxhet) { 	next; 	}
	
#warn "--------------\n$start[0]_$start[1]\n\n";
	my @rest;
	my $match=0;
	my $miss=0;
	my $anum=0;
	my $g1;
	my $g2;
	my $a1;
	my $a2;
	my $nalt=0;
	my $nref=0;
	my $het=0;
	my $nalleles=0;
	my $altpairs=0;
	my %seen={};
	for(my $p=0;$pp=$pairs[$p];$p++) {
		($r1,$r2)=split(":",$npairs[$p]);
		($i1,$i2)=split(":",$pp);
		($g1,@rest)=split(":",$lin[$i1]);
		($g2,@rest)=split(":",$lin[$i2]);
		if ($g1=~/\./ || $g2=~/\./) { 
			$miss++;
			$match++;
		}
		elsif ($g1 eq $g2) { 
			$match++;
			($a1,$a2)=split(/\D/,$g1);
			if (!$seen{$a1}) {
				$seen{$a1}=1;
				$nalleles++;
			}
			if (!$seen{$a2}) {
				$seen{$a2}=1;
				$nalleles++;
			}
			if ($a1=~/[123456789]/) { $nalt++; }
			else { $nref++;}
			if ($a2=~/[123456789]/) { 
				$nalt++;
				$altpairs++;
			}		
			else { $nref++;}
			if ($a1 ne $a2) {$het++;}
		}
	}
	next if ($match < ($nreps*$fmatch) );
	next if ( ($miss/$nreps) > $missing);
	$pass++;
	if ($nalt) { $numalt++;} else {next;}
	next if ($altpairs<$altPairs );
	$altpass++; 
	next if ( $het<$hetPairs );
	$hetpass++; 
	if ($nalleles>1) { 
		$poly++;
		$numout++;
		print join("\t",@start)."\t".join("\t",@lin)."\n";
	}
	elsif (!$polyonly){
		if ($nalt) {
			$numout++;
			print join("\t",@start)."\t".join("\t",@lin)."\n";
		}
	}
}
warn "

$total total SNPs
$pass pass hets and match filters
$numalt show non-reference alleles
$altpass have alterantive alleles in at least $altPairs replicate pair(s)
$hetpass have matching heterozygotes in at least $hetPairs replicate pair(s)
$poly polymorphic\n$numout written

";	
