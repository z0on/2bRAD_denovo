#!/usr/bin/perl

my $usage="

recalibrateSNPs.pl :

Non-parametric variant quality recalibration based on INFO fields 
in the set of variants that are reproducibly genotyped among replicates 
(output of replicatesMatch.pl)

This script is for de novo pipeline. Use recalibrateSNPs_gatk.pl for
reference-based pipeline.

Three parameters (INFO items) are assessed for each SNP: 
- total coverage depth (DP)
- strand bias (SB)
- allele bias (AB)
- SNP position withn tag (TP)

For AB and SB, low values are deemed bad and high ones are good;
for DP and TP, both  low and high values are considered bad, the best values
being around the median.

The script computes the product of quantiles for different INFO fields, determines 
the quantiles of the result within the 'true' set (JOINT quantiles),
and then computes the new quality scores in the main vcf file.

These scores are supposed to correspond to the probability (x100) that the 
SNPs comes from the same distribution as the 'true' SNPs. 

Output:

- Recalibrated VCF (printed to STDOUT) with QUAL field replaced by the new score 
  minus 0.1 (for easy filtering)
  
- a table (printed to STDERR) showing the \"gain\" at each recalibrated quality 
  score setting, which is excess variants removed when filtering at this quality score.
  For example, if 40% of all variants are removed at the quality score 15, the gain 
  is 40 - 15 = 35% (i.e., in addition to removing 15% of the 'true' variants, 
  additional 35% of the dataset, likely corresponding to wrong variants, is removed). 
  The optimal filtering score is the one giving maximum gain. Try different combinations 
  of the four possible filters to find the one maximizing the gain.

Arguments:

vcf=[file name]  : vcf file to be recalibrated

true=[file name] : vcf file that is subset of the above, with SNPs that are considered 
                   true because they show matching and polymorphic genotypes in replicates 
                   (replicatesMatch.pl polyonly=1). 
                   
           -nodp : do not use DP 
           -noab : do not use AB
           -nosb : do not use SB
           -notp : do not use TP

Example:
recalibrateSNPs.pl vcf=cdh_alltags.ul_Vhap_count10_ab10_sb10_clip0.vcf \
true=vqsr.vhap.vcf -nosb -notp >denovo.vhap.recal.vcf

Mikhail Matz, matz\@utexas.edu July 2013 - September 2014
 
";

my $vcf;
my $true;
my $nodp=0;
my $nosb=0;
my $noab=0;
my $notp=0;

if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1;}
else { die $usage; }
if ("@ARGV"=~/true=(\S+)/) { $true=$1;}
else { die $usage; }
if ("@ARGV"=~/-noab/) { $noab=1;}
if ("@ARGV"=~/-nosb/) { $nosb=1;}
if ("@ARGV"=~/-nodp/) { $nodp=1;}
if ("@ARGV"=~/-notp/) { $notp=1;}
if($nodp+$nosb+$noab+$notp==4) { die "no fields left to filter!\n";}

open TR, $true or die "cannot open the true set $true\n";

my @dp=();
my @ab=();
my @sb=();
my @tp=();
my @dline=();
my @iline=();

while (<TR>) {
	next if ($_=~/^#/);
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	foreach my $i (@iline){
		if ($i=~/AB=(\d+)/) { push @ab, $1;}
		elsif ($i=~/SB=(\d+)/) { push @sb, $1;}
		elsif ($i=~/DP=(\d+)/) { push @dp, $1;}
		elsif ($i=~/TP=(\d+)/) { push @tp, $1;}
	}
}
close TR;

@ab=sort {$a <=> $b} @ab;
@sb=sort {$a <=> $b} @sb;
@dp=sort {$a <=> $b} @dp;
@tp=sort {$a <=> $b} @tp;

my $total=$#dp+1;

my %abq;
my %sbq;
my %dpq;
my %dpq2;

$abq{1}=$ab[sprintf("%.0f",$total*0.01)];
$sbq{1}=$sb[sprintf("%.0f",$total*0.01)];
$abq{5}=$ab[sprintf("%.0f",$total*0.05)];
$sbq{5}=$sb[sprintf("%.0f",$total*0.05)];
$abq{10}=$ab[sprintf("%.0f",$total*0.1)];
$sbq{10}=$sb[sprintf("%.0f",$total*0.1)];
$abq{15}=$ab[sprintf("%.0f",$total*0.15)];
$sbq{15}=$sb[sprintf("%.0f",$total*0.15)];
$abq{20}=$ab[sprintf("%.0f",$total*0.2)];
$sbq{20}=$sb[sprintf("%.0f",$total*0.2)];
$abq{30}=$ab[sprintf("%.0f",$total*0.3)];
$sbq{30}=$sb[sprintf("%.0f",$total*0.3)];
$abq{40}=$ab[sprintf("%.0f",$total*0.4)];
$sbq{40}=$sb[sprintf("%.0f",$total*0.4)];
$abq{50}=$ab[sprintf("%.0f",$total*0.5)];
$sbq{50}=$sb[sprintf("%.0f",$total*0.5)];
$abq{60}=$ab[sprintf("%.0f",$total*0.6)];
$sbq{60}=$sb[sprintf("%.0f",$total*0.6)];
$abq{70}=$ab[sprintf("%.0f",$total*0.7)];
$sbq{70}=$sb[sprintf("%.0f",$total*0.7)];
$abq{80}=$ab[sprintf("%.0f",$total*0.8)];
$sbq{80}=$sb[sprintf("%.0f",$total*0.8)];
$abq{90}=$ab[sprintf("%.0f",$total*0.9)];
$sbq{90}=$sb[sprintf("%.0f",$total*0.9)];


$dpq{1}=$dp[sprintf("%.0f",$total*0.005)];
$dpq2{1}=$dp[sprintf("%.0f",$total*0.995)];
$dpq{5}=$dp[sprintf("%.0f",$total*0.025)];
$dpq2{5}=$dp[sprintf("%.0f",$total*0.975)];
$dpq{10}=$dp[sprintf("%.0f",$total*0.05)];
$dpq2{10}=$dp[sprintf("%.0f",$total*0.95)];
$dpq{15}=$dp[sprintf("%.0f",$total*0.075)];
$dpq2{15}=$dp[sprintf("%.0f",$total*0.925)];
$dpq{20}=$dp[sprintf("%.0f",$total*0.1)];
$dpq2{20}=$dp[sprintf("%.0f",$total*0.9)];
$dpq{30}=$dp[sprintf("%.0f",$total*0.15)];
$dpq2{30}=$dp[sprintf("%.0f",$total*0.85)];
$dpq{40}=$dp[sprintf("%.0f",$total*0.2)];
$dpq2{40}=$dp[sprintf("%.0f",$total*0.8)];
$dpq{50}=$dp[sprintf("%.0f",$total*0.25)];
$dpq2{50}=$dp[sprintf("%.0f",$total*0.75)];
$dpq{60}=$dp[sprintf("%.0f",$total*0.3)];
$dpq2{60}=$dp[sprintf("%.0f",$total*0.7)];
$dpq{70}=$dp[sprintf("%.0f",$total*0.35)];
$dpq2{70}=$dp[sprintf("%.0f",$total*0.65)];
$dpq{80}=$dp[sprintf("%.0f",$total*0.4)];
$dpq2{80}=$dp[sprintf("%.0f",$total*0.6)];
$dpq{90}=$dp[sprintf("%.0f",$total*0.45)];
$dpq2{90}=$dp[sprintf("%.0f",$total*0.55)];

$tpq{1}=$tp[sprintf("%.0f",$total*0.005)];
$tpq2{1}=$tp[sprintf("%.0f",$total*0.995)];
$tpq{5}=$tp[sprintf("%.0f",$total*0.025)];
$tpq2{5}=$tp[sprintf("%.0f",$total*0.975)];
$tpq{10}=$tp[sprintf("%.0f",$total*0.05)];
$tpq2{10}=$tp[sprintf("%.0f",$total*0.95)];
$tpq{15}=$tp[sprintf("%.0f",$total*0.075)];
$tpq2{15}=$tp[sprintf("%.0f",$total*0.925)];
$tpq{20}=$tp[sprintf("%.0f",$total*0.1)];
$tpq2{20}=$tp[sprintf("%.0f",$total*0.9)];
$tpq{30}=$tp[sprintf("%.0f",$total*0.15)];
$tpq2{30}=$tp[sprintf("%.0f",$total*0.85)];
$tpq{40}=$tp[sprintf("%.0f",$total*0.2)];
$tpq2{40}=$tp[sprintf("%.0f",$total*0.8)];
$tpq{50}=$tp[sprintf("%.0f",$total*0.25)];
$tpq2{50}=$tp[sprintf("%.0f",$total*0.75)];
$tpq{60}=$tp[sprintf("%.0f",$total*0.3)];
$tpq2{60}=$tp[sprintf("%.0f",$total*0.7)];
$tpq{70}=$tp[sprintf("%.0f",$total*0.35)];
$tpq2{70}=$tp[sprintf("%.0f",$total*0.65)];
$tpq{80}=$tp[sprintf("%.0f",$total*0.4)];
$tpq2{80}=$tp[sprintf("%.0f",$total*0.6)];
$tpq{90}=$tp[sprintf("%.0f",$total*0.45)];
$tpq2{90}=$tp[sprintf("%.0f",$total*0.55)];


print STDERR "\nAB quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %abq){
	print STDERR "$q\t$abq{$q}\n";
}

print STDERR "\nSB quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %sbq){
	print STDERR "$q\t$sbq{$q}\n";
}

print STDERR "\nDP quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %dpq){
	print STDERR "$q\t$dpq{$q}\t$dpq2{$q}\n";
}

print STDERR "\nTP quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %tpq){
	print STDERR "$q\t$tpq{$q}\t$tpq2{$q}\n";
}

my @joint=();
my %jq={};

open TR, $true;
while(<TR>){
	if ($_=~/^#/) { 
		next;
	}
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	my $abQ=0;
	my $sbQ=0;
	my $dpQ=0;
	my $tpQ=0;
	foreach my $i (@iline){
		if ($i=~/AB=(\d+)/) { 
			if ($noab){
				$abQ=100;
				next;
			}
			my $abt=$1;
			foreach my $abv (sort {$a <=> $b} keys %abq) {
				if ($abt<$abq{$abv}) { 
					$abQ=$abv;
					last;
				}
			}
			if (!$abQ){$abQ=100;}			
		}
		elsif ($i=~/SB=(\d+)/) { 
			if ($nosb){
				$sbQ=100;
				next;
			}
			my $sbt=$1;
			foreach my $sbv (sort {$a <=> $b} keys %sbq) {
				if ($sbt<$sbq{$sbv}) { 
					$sbQ=$sbv;
					last;
				}
			}			
			if (!$sbQ){$sbQ=100;}			
		}
		elsif ($i=~/DP=(\d+)/) { 
			if ($nodp){
				$dpQ=100;
				next;
			}
			my $dpt=$1;
			foreach my $dpv (sort {$a <=> $b} keys %dpq) {
				if ($dpt<$dpq{$dpv} || $dpt>$dpq2{$dpv}) { 
					$dpQ=$dpv;
					last;
				}
			}
			if (!$dpQ){
				$dpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
		}
		elsif ($i=~/TP=(\d+)/) { 
			if ($notp){
				$tpQ=100;
				next;
			}
			my $tpt=$1;
			foreach my $tpv (sort {$a <=> $b} keys %tpq) {
				if ($tpt<$tpq{$tpv} || $tpt>$tpq2{$tpv}) { 
					$tpQ=$tpv;
					last;
				}
			}
			if (!$tpQ){
				$tpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
		}
	}
	push @joint, sprintf("%.5f",($dpQ/100)*($tpQ/100)*($abQ/100)*($sbQ/100));
}
close TR;

@joint=sort {$a <=> $b} @joint;
$jq{1}=$joint[sprintf("%.0f",$total*0.01)];
$jq{5}=$joint[sprintf("%.0f",$total*0.05)];
$jq{10}=$joint[sprintf("%.0f",$total*0.1)];
$jq{15}=$joint[sprintf("%.0f",$total*0.15)];
$jq{20}=$joint[sprintf("%.0f",$total*0.2)];
$jq{30}=$joint[sprintf("%.0f",$total*0.3)];
$jq{40}=$joint[sprintf("%.0f",$total*0.4)];
$jq{50}=$joint[sprintf("%.0f",$total*0.5)];
$jq{60}=$joint[sprintf("%.0f",$total*0.6)];
$jq{70}=$joint[sprintf("%.0f",$total*0.7)];
$jq{80}=$joint[sprintf("%.0f",$total*0.8)];
$jq{90}=$joint[sprintf("%.0f",$total*0.9)];

print STDERR "\nJOINT quantiles:\n";
foreach my $q (sort {$a <=> $b} keys %jq){
	next if ($q=~/HASH/);
	print STDERR "$q\t$jq{$q}\n";
}

my $ones=0;
my $fives=0;
my $tens=0;
my $fifteens=0;
my $twenties=0;
my $thirties=0;
my $total2=0;

open VCF, $vcf or die "cannot open vcf $vcf\n";

while(<VCF>){
	if ($_=~/^#/) { 
		print $_;
		next;
	}
	$total2++;
	@dline=split("\t",$_);
	@iline=split(";",$dline[7]);
	my $abQ=0;
	my $sbQ=0;
	my $dpQ=0;
	my $jQ=0;
	foreach my $i (@iline){
		if ($i=~/AB=(\d+)/) { 
			if ($noab){
				$abQ=100;
				next;
			}
			my $abt=$1;
			foreach my $abv (sort {$a <=> $b} keys %abq) {
				if ($abt<$abq{$abv}) { 
					$abQ=$abv;
					last;
				}
			}
			if (!$abQ){$abQ=100;}			
		}
		elsif ($i=~/SB=(\d+)/) { 
			if ($nosb){
				$sbQ=100;
				next;
			}
			my $sbt=$1;
			foreach my $sbv (sort {$a <=> $b} keys %sbq) {
				if ($sbt<$sbq{$sbv}) { 
					$sbQ=$sbv;
					last;
				}
			}			
			if (!$sbQ){$sbQ=100;}			
		}
		elsif ($i=~/DP=(\d+)/) { 
			if ($nodp){
				$dpQ=100;
				next;
			}
			my $dpt=$1;
			foreach my $dpv (sort {$a <=> $b} keys %dpq) {
				if ($dpt<$dpq{$dpv} || $dpt>$dpq2{$dpv}) { 
					$dpQ=$dpv;
					last;
				}
			}
			if (!$dpQ){
				$dpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
		}
		elsif ($i=~/TP=(\d+)/) { 
			if ($notp){
				$tpQ=100;
				next;
			}
			my $tpt=$1;
			foreach my $tpv (sort {$a <=> $b} keys %tpq) {
				if ($tpt<$tpq{$tpv} || $tpt>$tpq2{$tpv}) { 
					$tpQ=$tpv;
					last;
				}
			}
			if (!$tpQ){
				$tpQ=100;
			}
#			if ($dpQ>50) { $dpQ=100-$dpQ;}						
		}
	}
	my $jt=sprintf("%.5f",($dpQ/100)*($tpQ/100)*($abQ/100)*($sbQ/100));
#warn "\t$jt\n";
	foreach my $jv (sort {$a <=> $b} keys %jq) {
		if ($jt<$jq{$jv}) { 
			$jQ=$jv;
			last;
		}
	}			
	if (!$jQ){$jQ=100;}
	if ($jQ==1) { $ones++;}
	elsif ($jQ==5) { $fives++;}
	elsif ($jQ==10) { $tens++;}			
	elsif ($jQ==15) { $fifteens++;}			
	elsif ($jQ==20) { $twenties++;}			
	elsif ($jQ==30) { $thirties++;}	
	$jQ=$jQ-0.1;		
	$dline[5]=$jQ;
	push @iline,"SBQ=".$sbQ;
	push @iline,"ABQ=".$abQ;
	push @iline,"DPQ=".$dpQ;
	push @iline,"TPQ=".$tpQ;
	$dline[7]=join(";",@iline);	
	print join("\t",@dline);
}

my $oness=sprintf("%.2f",100*$ones/$total2);
my $fivess=sprintf("%.2f",100*($ones+$fives)/$total2);
my $tenss=sprintf("%.2f",100*($ones+$fives+$tens)/$total2);
my $fifteenss=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens)/$total2);
my $twentiess=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens+$twenties)/$total2);
my $thirtiess=sprintf("%.2f",100*($ones+$fives+$tens+$fifteens+$twenties+$thirties)/$total2);

print STDERR "
------------------------
$oness%\tat qual <1 (",sprintf("%.2f",$oness-1),"% gain)
$fivess%\tat qual <5 (",sprintf("%.2f",$fivess-5),"% gain)
$tenss%\tat qual <10 (",sprintf("%.2f",$tenss-10),"% gain)
$fifteenss%\tat qual <15 (",sprintf("%.2f",$fifteenss-15),"% gain)
$twentiess%\tat qual <20 (",sprintf("%.2f",$twentiess-20),"% gain)
$thirtiess%\tat qual <30 (",sprintf("%.2f",$thirtiess-30),"% gain)
------------------------

";


