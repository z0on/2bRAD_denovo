#!/usr/bin/perl

my $usage="
thinner.pl (v.4.2) : 

Leaves only one SNP per specified nucleotide distance interval, chosen by its 
allele frequency, sequencing depth, combination of those, or randomly.

The minimal distance between selected SNPs is half the specified interval.

The thinning process uses random starts for interval specification, potentially
resulting in different SNPs selected. 

Arguments:

               vcf=[file name] : vcf file name
            
            interval=[integer] : interval length, default 40 (for 2bRAD tags)
            
criterion=[maxAF|maxDP-random|maxDP-maxAF|random] : SNP choosing criterion.
                             maxAF - maximum minor allele frequency;
                      maxDP-random - randomly selected from those with maximum 
                                     sequencing depth (for dadi);
                       maxDP-maxAF - maximum allele frequency among those with
                                     maximum sequencing depth (default);
                            random - random selection.

Output:  thinned VCF, printed to STDOUT

Example: thinner.pl vcf=denovo.recal.vcf > thinDenov.vcf

NOTES: - Do not thin variants if you want to calculate pi or Tajima's D.

Mikhail Matz, matz\@utexas.edu

";

my $vcf;
if (" @ARGV "=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage; }
my $inter=80;
if (" @ARGV "=~/interval=(\d+)/) { $inter=$1*2; } 
my $criterion="maxDP-maxAF";
if (" @ARGV "=~/criterion=maxAF/) { $criterion="maxAF"; } 
elsif (" @ARGV "=~/criterion=random/) { $criterion="random"; } 
elsif (" @ARGV "=~/criterion=maxDP-random/) { $criterion="maxDP-random"; } 

#warn "criterion:$criterion\n";

open VCF, $vcf or die "cannot open vcf file $vcf\n";

my @dats;
my $info;
my $chrom;
my $pos;
my $pos0=-sprintf("%.0f",rand $inter/2);
my $poslast=-$inter;
#warn "pos0:$pos0;poslast:$poslast\n";
my $af;
my $dp;
my @afs=();
my @dps=();
my $maxaf=0;
my $maxdp=0;
my @snps=();
my $chrom0="";
my @maxdps=();
my $skipped=0;
my $selected=0;
my $total=0;

while (<VCF>) {
	if ($_=~/^#/) { print and next;}
	chomp;
	$total++;
	@dats=split(/\t/,$_);
	$chrom=$dats[0];
	if (!$chrom0) { $chrom0=$chrom;}
	$pos=$dats[1];
	if ( $pos-$poslast < $inter/2 and ($chrom eq $chrom0)) {
#warn "\t\t\tskipping $pos\n" ;
		$skipped++;
		next;
	}
	$info=$dats[7];
	if ($info=~/AF=(\S*?)[:;]/) { $af=$1; } else { warn "no AF:\n@dats\n" and next;}
	if ($info=~/DP=(\d*?);/) { $dp=$1; } else { warn "no DP:\n@dats\n" and next;}
	if ($af>0.5) { $af=1-$af; }
	if (($pos-$pos0 < $inter) and ($chrom eq $chrom0)){
#warn "\t pos:$pos - interval:",$pos0,"-",$pos0+$inter," - poslast:$poslast\n";
		push @snps,join("\t",@dats);
		push @afs,$af;
		push @dps,$dp;
		if ($af >$afs[$maxaf]) { $maxaf=$#afs; }
		if ($dp > $dps[$maxdp]) { 
			@maxdps=();
			$maxdp=$#dps;
			push @maxdps, $maxdp; 
			}
		elsif ($dp==$dps[$maxdp]) { push @maxdps, $#dps; }
	}
	else {
		if ($snps[0]){
			my $datum;
#warn "$chrom0: $#snps | @afs | @dps :  maxaf:$afs[$maxaf]; maxdp:$dps[$maxdp] (@maxdps) \n";
			if ($criterion eq "random") { 
				$datum=$snps[rand @snps];
			}
			elsif ($criterion eq "maxAF") { 
				$datum=$snps[$maxaf] ;
			}
			elsif ($criterion eq "maxDP-random"){ 
				$datum=$snps[$maxdps[rand @maxdps]] ; 
			}
			elsif ($criterion eq "maxDP-maxAF"){
				$maxaf=0;
				foreach my $md (@maxdps) { 
					if ($afs[$md]>$maxaf) {
						$maxaf=$afs[$md];
						$maxdp=$md;
					}
				} 	
				$datum=$snps[$maxdp] ; 
#warn "chose $maxdp\n";
			}
			my @datss=split(/\t/,$datum);
			$poslast=$datss[1];
#warn "\t\tSel:$poslast\n";
			print $datum,"\n" ;
			$selected++; 				 
		}
		$maxaf=0;
		$maxdp=0;
		@afs=();
		@dps=();
		@snps=();
		if ( $pos-$poslast < $inter/2 and ($chrom eq $chrom0)) {
#warn "\t\t\tskipping $pos\n" ;
			$skipped++;
			next;
		}
		if ($chrom ne $chrom0) {
			$chrom0=$chrom;
#warn "\n-------------\nnew chrom: $chrom\n";
			$poslast=-$inter;
			$pos0=-sprintf("%.0f",rand $inter/2);
		}	
		while ($pos0+$inter<$pos) { $pos0+=$inter;}
#warn "\t Pos:$pos - interval:",$pos0,"-",$pos0+$inter," - poslast:$poslast\n";
		push @afs,$af;
		push @dps,$dp;
		@maxdps=();
		push @maxdps, 0; 
		push @snps,join("\t",@dats);
	}
}

if ($snps[0]){
	my $datum;
#warn "$chrom0: $#snps | @afs | @dps :  maxaf:$afs[$maxaf]; maxdp:$dps[$maxdp] (@maxdps) \n";
	if ($criterion eq "random") { 
		$datum=$snps[rand @snps];
	}
	elsif ($criterion eq "maxAF") { 
		$datum=$snps[$maxaf] ;
	}
	elsif ($criterion eq "maxDP-random"){ 
		$datum=$snps[$maxdps[rand @maxdps]] ; 
	}
	elsif ($criterion eq "maxDP-maxAF"){
		$maxaf=0;
		foreach my $md (@maxdps) { 
			if ($afs[$md]>$maxaf) {
				$maxaf=$afs[$md];
				$maxdp=$md;
			}
		} 	
		$datum=$snps[$maxdp] ; 
#warn "chose $maxdp\n";
	}
	my @datss=split(/\t/,$datum);
	$poslast=$datss[1];
#warn "\t\tSel:$poslast\n";
	print $datum,"\n" ;
	$selected++; 				 
}

warn "
$total total loci
$skipped loci skipped because they were closer than ",$inter/2,"
$selected loci selected

";