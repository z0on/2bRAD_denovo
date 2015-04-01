#!/usr/bin/perl

my $usage="

repMatchStats.pl :

Summarizes genotypic match between replicates in a vcf file.

Arguments:

       vcf=[file name] : input vcf file
replicates=[file name] : a two column tab-delimited table listing 
                         pairs of samples that are replicates

Output:
A table (printed to STDOUT) of the following form:

pair      gtyped  match      [ 00  01  11 ] HetMatch HomoHetMism HetNoCall HetsDiscRate
K210:K212  7328  7169(97.8%) [ 78% 17% 5% ]  1200         139	     5         0.94	
K212:K213  7369  7117(96.6%) [ 78% 17% 5% ]  1202         179        4         0.93	

The first four columns are self-evident;
the last four columns show how good is the match between heterozygote calls 
(HetNoCall being a heterozygote in one and missing data in another replicate).
The last column is the most important: it is the heterozygote discovery rate,
the fraction of all heterozygotes discovered in each replicate.

Example: repMatchStats.pl vcf=denovo.filt.recode.vcf replicates=clonepairs.tab

Mikhail V. Matz, matz\@utexas.edu, July 2013 

";

my $vcf;
my $reps;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/replicates=(\S+)/) { $reps=$1;}
else { die $usage; }

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @samples;
my @pairs=();
my %indr={};
my @npairs=();
my $r1;
my $r2;
my $s;

my $nreps=0;

my @gtyped=();
my @match=();
my @homomatch=();
my @refhomomatch=();
my @althomomatch=();
my @heteromatch=();
my @mismatch=();
my @homohetmismatch=();
my @allhetmismatch=();
my @nocallhomo;
my @nocallhet;
my @misDP=();
my $dpslot=0;

while (<VCF>) {
	if ($_=~/^#/) {
		if ($_=~/contig/) { next;}
		elsif ($_=~/^#CHROM/){
#			print $_;
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
#warn "$nreps replicate pairs\n";
		}
#		else { print $_;}
		next;
	}
	chop;
	$total++;
	my @lin=split("\t",$_);
	if (!$dpslot){
		my $info=$lin[8];
		my @ifields=split(":",$info);
		for(my $ifi=0;$ifi<=$#ifields;$ifi++){
			if ($ifields[$ifi] eq "DP"){
				$dpslot=$ifi;
				last;
			}
		}
#print "@ifields    => DP:$dpslot\n";
	}
	splice(@lin,0,9);
#warn "--------------\n$start[0]_$start[1]\n\n";
	my @rest1;
	my @rest2;
	my $g1;
	my $g2;
	my $a1;
	my $a2;
	my $b1;
	my $b2;
	for(my $p=0;$pp=$pairs[$p];$p++) {
		($r1,$r2)=split(":",$npairs[$p]);
		($i1,$i2)=split(":",$pp);
		($g1,@rest1)=split(":",$lin[$i1]);
		($g2,@rest2)=split(":",$lin[$i2]);
#warn "$r1:$i1\t$g1\t$r2:$i2\t$g2\n";
#		next if ($g1=~/\./ || $g2=~/\./);
		($a1,$a2)=split(/[\/\|]/,$g1);
		($b1,$b2)=split(/[\/\|]/,$g2);
		if ($g1 eq $g2) { 
			next if ($a1=~/\./ && $a2=~/\./);
			$gtyped[$p]++;
			$match[$p]++;
			if ($a1==$a2) { 
				$homomatch[$p]++;
				if ($a1==0) { $refhomomatch[$p]++;}
				else { $althomomatch[$p]++;}
			}
			else { $heteromatch[$p]++;}
		}
	 	else {
	 		if ($a1=~/\./ && $a2=~/\./) {
	 			$nocall[$p]++;
	 			if ($b1==$b2) { $nocallhomo[$p]++;}
	 			else { $nocallhet[$p]++;}
	 		} 
	 		elsif ($b1=~/\./ && $b2=~/\./) {
	 			$gtyped[$p]++;
	 			$nocall[$p]++;
	 			if ($a1==$a2 ) { $nocallhomo[$p]++;}
	 			else { $nocallhet[$p]++;}
	 		} 
			else {
				$gtyped[$p]++;
				$mismatch[$p]++;
				if ($a1==$a2 ) {
					if ($b1==$b2) { 
						$allhomomismatch[$p]++;
					}
					else {
						$homohetmismatch[$p]++;
						push @misDP, $rest1[$dpslot-1];
					}
				}
				else {
					if ($b1==$b2 ) { 
						$homohetmismatch[$p]++;
						push @misDP, $rest2[$dpslot-1];
					}
					else { $allhetmismatch[$p]++; }
				}		
			}
		}
	}
}

@misDP=sort {$a <=> $b} @misDP;

my $mdp25=$misDP[sprintf("%.0F",($#misDP+1)*0.25)];
my $mdp50=$misDP[sprintf("%.0F",($#misDP+1)*0.5)];
my $mdp75=$misDP[sprintf("%.0F",($#misDP+1)*0.75)];
 
print "pair\tgtyped\tmatch\t[ 00\t01\t11 ]\tHetMatch\tHomoHetMismatch\tHetNoCall\tHetsDiscoveryRate\n";
for(my $p=0;$pp=$npairs[$p];$p++) {
	print "$pp\t$gtyped[$p]\t$match[$p](",sprintf("%.1f",100*$match[$p]/$gtyped[$p]),"%)\t [";
	print sprintf("%.0f",100*$refhomomatch[$p]/$match[$p]),"%\t";
	print sprintf("%.0f",100*$heteromatch[$p]/$match[$p]),"%\t";
	print sprintf("%.0f",100*$althomomatch[$p]/$match[$p]),"% ]\t";
	print $heteromatch[$p],"\t";
	print $homohetmismatch[$p],"\t";
	print $nocallhet[$p],"\t";
#	print $nocallhomo[$p],"\t";
#	print $allhetmismatch[$p],"\t";
#	print $allhomomismatch[$p],"\t";
	my $matchfrac=2*$heteromatch[$p]/($nocallhet[$p]+2*$heteromatch[$p]+$homohetmismatch[$p]);
	print sprintf("%.2f",$matchfrac),"\t\n";

#	print sprintf("%.1f",100*($refhomohetmismatch[$p]+$althomohetmismatch[$p])/($refhomohetmismatch[$p]+$althomohetmismatch[$p]+$heteromatch[$p])),"%\t";
#	print "$mismatch[$p](",sprintf("%.1f",100*$mismatch[$p]/$gtyped[$p]),"%)\t[ ";
#	print sprintf("%.1f",100*$allhetmismatch[$p]/$mismatch[$p]),"%\t";
#	print sprintf("%.1f",100*$homohetmismatch[$p]/$mismatch[$p]),"%\t[ ";
#	print sprintf("%.1f",100*$refhomohetmismatch[$p]/$homohetmismatch[$p]),"%\t";
#	print sprintf("%.1f",100*$althomohetmismatch[$p]/$homohetmismatch[$p]),"% ]]\n";
}

print "
------------------------
hets called homos depth: 
lower 25%	$mdp25
median		$mdp50
upper 75%	$mdp75

";