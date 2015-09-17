#!/usr/bin/perl

my $usage="

thinner.pl (v.3) : 

Leaves only one SNP per specified nucleotide distance interval, chosen by its 
allele frequency, sequencing depth, combination of those, or randomly.

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

Mikhail Matz, matz@utexas.edu

";

my $vcf;
if (" @ARGV "=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage; }
my $inter=40;
if (" @ARGV "=~/interval=(\d+)/) { $inter=$1; } 
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
my $pos0=-1000000;
my $af;
my $dp;
my @afs=();
my @dps=();
my $maxaf=0;
my $maxdp=0;
my @snps=();
my $chrom0="";
my @maxdps=();

while (<VCF>) {
	if ($_=~/^#/) { print and next;}
	chomp;
	@dats=split(/\t/,$_);
	$chrom=$dats[0];
	$pos=$dats[1];
	$info=$dats[7];
	if ($info=~/AF=(\S*?);/) { $af=$1; } else { warn "no AF:\n@dats\n" and next;}
	if ($info=~/DP=(\d*?);/) { $dp=$1; } else { warn "no DP:\n@dats\n" and next;}
	if ($af>0.5) { $af=1-$af; }
	if ($pos-$pos0 <=$inter and ($chrom eq $chrom0)){
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
#warn "$pos0: $#snps | @afs | @dps :  maxaf:$afs[$maxaf]; maxdp:$dps[$maxdp] (@maxdps) \n";
			if ($criterion eq "random") { print $snps[rand @snps],"\n" ; }
			elsif ($criterion eq "maxAF") { print $snps[$maxaf] ,"\n"  ; }
			elsif ($criterion eq "maxDP-random"){ print $snps[$maxdps[rand @maxdps]] ,"\n" ; }
			elsif ($criterion eq "maxDP-maxAF"){
				$maxaf=0;
				foreach my $md (@maxdps) { 
					if ($afs[$md]>$maxaf) {
						$maxaf=$afs[$md];
						$maxdp=$md;
					}
				} 	
				print $snps[$maxdp] ,"\n" ; 
#warn "chose $maxdp\n";
			}
		}
		$pos0=$pos;
		$chrom0=$chrom;
		$maxaf=0;
		$maxdp=0;
		@afs=();
		@dps=();
		@snps=();
		push @afs,$af;
		push @dps,$dp;
		@maxdps=();
		push @maxdps, 0; 
		push @snps,join("\t",@dats);
	}
}

if ($snps[0]){
#warn "$pos0: $#snps | @afs | @dps :  maxaf:$afs[$maxaf]; maxdp:$dps[$maxdp] (@maxdps) \n";
	if ($criterion eq "random") { print $snps[rand @snps],"\n" ; }
	elsif ($criterion eq "maxAF") { print $snps[$maxaf] ,"\n"  ; }
	elsif ($criterion eq "maxDP-random"){ print $snps[$maxdps[rand @maxdps]] ,"\n" ; }
	elsif ($criterion eq "maxDP-maxAF"){
		$maxaf=0;
		foreach my $md (@maxdps) { 
			if ($afs[$md]>$maxaf) {
				$maxaf=$afs[$md];
				$maxdp=$md;
			}
		} 	
		print $snps[$maxdp] ,"\n" ; 
#warn "chose $maxdp\n";
	}
}
