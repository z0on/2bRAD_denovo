#!/usr/bin/perl

my $usage="

replicatesMatch.pl : (version 0.3, April 14, 2015)

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

allAlts=[1|0]    output all SNPs showing non-reference alleles (not necessarily
                 polymorphic among genotyped samples). Default 0

max.het=[float]  maximum fraction of heterozygotes (guards against lumped paralogous 
                 loci). Default 0.75.

Example: 
replicatesMatch.pl vcf=cdh_alltags.ul_Variants_count10_ab10_sb10_clip0.vcf \
replicates=clonepairs.tab polyonly=1 > vqsr.denovo.vcf

";

my $vcf;
my $reps;
my $missing=0.25;
my $fmatch=1;
my $altonly=0;
my $maxhet=0.75;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/replicates=(\S+)/) { $reps=$1;}
else { die $usage; }
if ("@ARGV"=~/missing=(\S+)/) { $missing=$1;}
if ("@ARGV"=~/matching=(\S+)/) { $fmatch=$1;}
if ("@ARGV"=~/allAlts=1/ ) { $altonly=1;}
if ("@ARGV"=~/max.het=(\S+)/) { $maxhet=$1;}

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @samples;
my @pairs=();
my %indr={};
my @npairs=();
my $r1;
my $r2;
my $nreps=0;
my $pass=0;
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
#warn "--------------\n$start[0]_$start[1]\n\n";
	my @rest;
	my $match=0;
	my $miss=0;
	my %seen={};
	my $anum=0;
	my $g1;
	my $g2;
	my $a1;
	my $a2;
	my $nalt=0;
	my $nref=0;
	my $het=0;
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
			if ($a1=~/[123456789]/) { $nalt++;}
			else { $nref++;}
			if ($a2=~/[123456789]/) { $nalt++;}			
			else { $nref++;}
			if ($a1 ne $a2) {$het++;}
		}
	}
	next if ($match < ($nreps*$fmatch) );
	next if ( ($miss/$nreps) > $missing);
	if ($nalt) { $numalt++;}
	if ($nref and $nalt) { 
		$poly++;
		if ($het/$nreps<$maxhet) { 	
			$numout++;
			print join("\t",@start)."\t".join("\t",@lin)."\n";

#warn "--------------\n$start[0]_$start[1]\n";
#for(my $p=0;$pp=$pairs[$p];$p++) {
#	($r1,$r2)=split(":",$npairs[$p]);
#	($i1,$i2)=split(":",$pp);
#	($g1,@rest)=split(":",$lin[$i1]);
#	($g2,@rest)=split(":",$lin[$i2]);
#warn "$r1:$i1\t$g1\t$r2:$i2\t$g2\n";
#}
#warn "het:$het nalt:$nalt nref:$nref miss:$miss match:$match\n";

		}
	}
	elsif ($altonly){
		if ($nalt and $het/$nreps<$maxhet) {
			$numout++;
			print join("\t",@start)."\t".join("\t",@lin)."\n";
		}
	}
	$pass++;
}
warn "$total total SNPs\n$pass passed\n$numalt alt alleles\n$poly polymorphic\n$numout written\n\n";	