#!/usr/bin/perl

my $usage="

haplocall_denovo.pl :

Calls genotypes (as population-wide RAD-tag haplotypes) in de novo 2bRAD data.

Outputs two files in VCF v.4 format: 
    \"Vhap\": whole RAD tags recoded as alternative alleles (for coalescent analysis
              and methods that take advantage of multi-allelic markers);
\"Variants\": RAD tags broken into individual SNPs, within-tag phase recorded 
              (for methods requiring biallelic markers).

The QUAL field in VCF files is allele bias, the average fraction of minor allele 
counts in a heterozygote. Ideally this should be close to 50, more if the allele is
found in homozygotes as well.  

Recorded variant info fields:
SB : Strand bias (ratio x100)
AB : Minor allele bias in a heterozygote (ratio x100)
TP : SNP position in the tag
NS : Number of samples with data
DP : Combined read depth
AF : Allele frequency

Individual genotype fields:
GT : Genotype
GQ : Phred-scaled genotype quality
PS : Phase set (for multiple SNPs within a RAD tag)
AD : Allele depth
DP : Read depth

Homozygote quality is calculated as phred-scaled
p(error) = 0.5^(n-1), where n is the number of reads;
Heterozygote quality is the same plus 20 on phred-scale 
(we are very sure about our heterozygotes, see allele filters below)

The ends of tags with more than one SNP in the terminal 3 bases
are clipped off since such SNPs are most likely due to indel near the
end of the tag.

Required arguments and defaults:

arg1              input file (produced by uniq2loci.pl)

Optional arguments:

clip=0            number of bases to clip off the ends of reads

----- Allele filters -----

aobs=2            min number of times an allele must be observed across all samples. 

strbias=10        strand bias cutoff: minimal ratio between reverse and 
                  direct reads, x100. 0 would mean no filtering.
                  
abias=10          allele bias cutoff: minimal fraction of reads per locus corresponding 
                  to the candidate allele, averaged across individuals
                  containing the allele, x100. 

ind=1             minimum number of individuals in which a candidate new allele
                  has to be found.

----- Locus filters -----

found=0.25        fraction of individuals in which a locus has to be genotyped. 
                  Filter for this parameter later with vcftools.

hetero=0.5        maximum allowed fraction of heterozygotes. Guards against lumped 
                  paralogous loci.

mono=toss         whether to keep non-polymorphic loci. Specify mono=keep to keep them
                  (if coalescent analysis with MIGRATE is planned downstream).

----- Genotype filters -----

mindp=3           minimum depth to call a homozygote. 


Mikhail Matz, July 2013 - April 2015
matz\@utexas.edu

Example:  haplocall_denovo.pl cdh_alltags.ul

";

my $gmap=shift or die $usage;

my $cut=5;
my $strb=10;
my $indfrac=0.25;
my $frachetero=0.5; 
my $monokeep=0;
my $numind=2;
my $minQ=10;
my $clip=0;
my $mindp=3;
if ("@ARGV"=~/strbias=(\d+)/) { $strb=$1;}
if ("@ARGV"=~/aobs=(\d+)/) { $cut=$1;}
if ("@ARGV"=~/found=(\S+)/) { $indfrac=$1;}
if ("@ARGV"=~/mono=keep/) { $monokeep=1;}
if ("@ARGV"=~/hetero=(\S+)/) { $frachetero=$1;}
if ("@ARGV"=~/ind=(\d+)/) { $numind=$1;}
if ("@ARGV"=~/abias=(\d+)/) { $minQ=$1;}
if ("@ARGV"=~/clip=(\d+)/) { $clip=$1;}
if ("@ARGV"=~/mindp=(\d+)/) { $mindp=$1;}

print "\nInput parameters (parameter names):
clipping $clip bases from read ends (clip)
must see an allele at least $cut times (aobs)
strand bias cutoff (strbias): $strb
allele bias cutoff (abias): $minQ
keep loci genotyped in at least $indfrac fraction of all samples (found)
must see an allele in at least $numind individual(s) (ind)
keep monomorphic tags? (mono): $monokeep 
maximum acceptable fracton of heterozygotes at a locus (hetero): $frachetero
minimum depth to call a homozygote (mindp): $mindp
\n";

my %seen={};
my $read;
my $ref;
my @rest=();
my %map={};
my $strand="";
my $seen2={};
my %locus={};
my %mapind={};
my %ot={};
my $now;
my @samples;
my $ct;
my $rc;
my @scounts;
my %sb={};
my %seen0={};

print "\tfile $gmap\n";
$now=localtime;
warn "reading start: $now\n";
open INP, $gmap or die "cannot open input file $gmap\n";
#	my %ot={};
while (<INP>) {
	chop;
	if ($_=~/revcom/) {
		@samples=split("\t",$_);
		splice(@samples, 0, 5);
		next;
	}
	($ref,$read,$seq,$ct,$rc,@scounts)=split('\t',$_);
	if ($clip) { $seq=substr $seq, $clip,(-1)*$clip; }
#print "read|$read\t$seq\tref|$ref\n";
	$seen{$seq}=$ct;
	$seen0{$seq}++;
	for($i=0;$sa=$samples[$i];$i++) {
		${$seen2{$sa}}{$seq}+=$scounts[$i];
	}
	if (!${$mapind{$ref}}{$seq}) {
		push @{$map{$ref}},$seq;
		${$mapind{$ref}}{$seq}=1;
	}
	my $d=$ct-$rc;
	if ($d>=$rc) { $sb{$seq}+=$rc/$d;}
	else {  $sb{$seq}+=$d/$rc;}
}
close INP;

foreach $seq (keys %seen0){
	next if ($seq=~/HASH/);
	$sb{$seq}=$sb{$seq}/$seen0{$seq};
	$sb{$seq}=sprintf("%.0f",100*$sb{$seq});
	undef $seen0{$seq};
}

#$now=localtime;
#warn "ref sorting start: $now\n";

my @refs=keys(%map);

$now=localtime;
warn "\nallele filters start: $now\n";

my $nsamples=$#samples+1;
my @goodrefs;
my %goodmap;
my %ni={};	

my $allcounts=0;
my $cutcounts=0;
my $indcounts=0;
my $strbcounts=0;
# excluding loci and alleles seen too few times
my $allloci=0;
my $cutloci=0;
foreach $r (@refs){
		$allloci++;
		next if ($r=~/HASH/);
		my $good=0;
		foreach $ss (@{$map{$r}}) {
			$allcounts++;
			if ($seen{$ss}<$cut) { 
				undef $seen{$ss};
				undef $sb{$ss};
				next;
			}
			$cutcounts++;
			foreach $sa (@samples) {
				if (${$seen2{$sa}}{$ss}) { 
					$ni{$ss}++;
				}
			}
			if ($ni{$ss}<$numind) { next; }
			$indcounts++;
			if ($sb{$ss}<$strb) { next;}
			$strbcounts++;
			foreach $sa (@samples) {
				if (${$seen2{$sa}}{$ss}) { 
						${$locus{$sa}}{$r}=1;
				}
			}
			push @{$goodmap{$r}},$ss;
			$good=1;
		}
		if ($good) {
			push @goodrefs, $r;
		}
}

print "\n-----------------\n
Allele filters:

$allcounts	raw alleles
$cutcounts	with $cut or more reads
$indcounts	in $numind or more samples
$strbcounts	pass strand bias cutoff $strb";

# measuring alleles' quality - 
# average fraction of tags in an individual containing the tag, *100
# excluding alleles with quality lower than $minQ 

my %goodmapQ={};
my $newcount=0;
foreach $r (@goodrefs){
	next if ($r=~/HASH/);
	my %aqual={};
	my %numsa={};
	foreach $sa (@samples) {
		my $tot=0;
		foreach $ss (@{$goodmap{$r}}) {
			next unless (${$seen2{$sa}}{$ss});
			$tot+=${$seen2{$sa}}{$ss};
		}
		my $QQ=0;
		foreach $ss (@{$goodmap{$r}}) {
			next unless (${$seen2{$sa}}{$ss});
			$numsa{$ss}++;
			$QQ=${$seen2{$sa}}{$ss}/$tot;
#			if ($QQ>0.5) { $QQ=1-$QQ;}
			$QQ=sprintf("%.2f",$QQ);
			$aqual{$ss}+=$QQ*100;
		}
	}
	foreach $ss (@{$goodmap{$r}}) {
		$aqual{$ss}=sprintf("%.0f",$aqual{$ss}/$numsa{$ss});
		next unless ($aqual{$ss}>=$minQ);
		push @{$goodmapQ{$r}},$ss;
		$newcount++;
	}
	$cutloci++; 	
}
print "
$newcount	pass allele bias cutoff $minQ
";

$now=localtime;
warn "\nlocus filters start: $now\n";

# excluding loci which have more than 2 alleles in any individual
my @goodrefs2=();
my $biallelic=0;
my $hetpass=0;
my $allloci=0;
foreach $r (@goodrefs){
	next if ($r=~/HASH/);
	my $nalleles=0;
	my $homo=0;
	my $genotyped=0;
	foreach $sa (@samples){
    		next if ($sa=~/HASH/);
		if (${$locus{$sa}}{$r}) {
		    $nalleles=0;	
		    $genotyped++;
			foreach $ss (@{$goodmapQ{$r}}) {
				if (${$seen2{$sa}}{$ss}) {
					$nalleles++;
				}
			}
			if ($nalleles>2) { last;}
			if ($nalleles==1) {$homo++;}
		}
	}
	$allloci++;
	next if ($homo<$genotyped*(1-$frachetero)); 
	$hetpass++;
	if ($nalleles>2) { next;}
	$biallelic++;
	push @goodrefs2, $r;
}

print "\n--------------------\n
Locus filters:

$allloci	total
$cutloci	remain after applying allele filters
$hetpass	have less than $frachetero fraction of heterozygotes
$biallelic 	with not more than 2 alleles per individual
";

# excluding loci genotyped in too few individuals
my @goodrefs3=();
foreach $r (@goodrefs2){
		my $gtyped=0;
		foreach $sa (@samples) {
				if (${$locus{$sa}}{$r}==1){ $gtyped++;}
		}
		if ($gtyped/$nsamples>=$indfrac) {
				push @goodrefs3, $r;
		}
}

print $#goodrefs3+1, "	genotyped in ",$indfrac*100,"% of samples\n";

my @goodrefs=();
if ($monokeep==1) { 
	@goodrefs=@goodrefs3;
	goto LOC;
}

# excluding monomorphic loci
foreach $r (@goodrefs3){
	next if ($r=~/HASH/);
	if ($#{$goodmapQ{$r}}<1) { next;}
	push @goodrefs, $r;
}

print $#goodrefs+1, "	polymorphic\n\n";

LOC:

$now=localtime;
warn "VCF writing start: $now\n";

# writing VCF format
$tabname=$gmap."_Variants_count".$cut."_ab".$minQ."_sb".$strb."_clip".$clip.".vcf";
open TAB, ">$tabname";
print {TAB} "##fileformat=VCFv4.1
##INFO=<ID=SB,Number=1,Type=Integer,Description=\"strand bias (ratio x100)\">
##INFO=<ID=AB,Number=1,Type=Integer,Description=\"minor allele bias (ratio x100)\">
##INFO=<ID=TP,Number=1,Type=Integer,Description=\"SNP position in the tag\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"combined read depth\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled genotype quality\">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele depth\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach $sa (@samples) {
	print {TAB} "\t$sa";
}
print {TAB} "\n";

my $chrom;
my $pos;
my $phase=0;
my $isphase=0;
my $tagprint=0;

my %goodsite={};
my %snps={};
my %trimtag={};

foreach $r (@goodrefs){
#print "ref $r\n";
	next if ($r=~/HASH/);
	$chrom=$r;
	my $maj="";
	my $maxcount=0;
	my $sumcount=0;
	foreach $ss (@{$goodmapQ{$r}}) {
		if ($seen{$ss}>$maxcount){
			$maxcount=$seen{$ss};
			$maj=$ss;
		}
		$sumcount+=$seen{$ss};
	}
	my @ordmap=();	
	push @ordmap, $maj;
	my $sbr=0;
	foreach $ss (@{$goodmapQ{$r}}) {
		$sbr+=$sb{$ss};
		next if ($ss eq $maj);
		push @ordmap, $ss;
	}
	$sbr=sprintf("%.0f",$sbr/($#{$goodmapQ{$r}}+1));
	my %gt={};
	my %gta={};
	foreach $s (@samples) {
		for ($a=0;$ssa=$ordmap[$a];$a++) {
			if(${$seen2{$s}}{$ssa}) {
				push @{$gt{$s}}, $a;
				push @{$gta{$s}},$ssa;
			}
		}
	}
#print "$maj major, $maxcount/$sumcount\n";	
	my @refbases=split("",$maj);
	my %snp={}; 
	my @ksnp;
	for ($b=0;$refbases[$b];$b++){
		foreach $ss (@ordmap) {
			my @abases=split ("",$ss);
			if ($refbases[$b] ne $abases[$b]) {
				push @ksnp,$b unless (" @ksnp "=~/ $b /);
#print "$ss : snp at ",$b+1,"\n";
				if (!$snp{$b}) { $snp{$b}=$refbases[$b];}
				$snp{$b}.=$abases[$b] unless ($snp{$b}=~/$abases[$b]/);
			}
		}
	}
#	my @ksnp=sort {$a <=> $b} keys %snp;	
	next if ("@ksnp"=~/12|13|22|23/);
	my $trim=0;
	if ($#ksnp>0) { 
#print STDERR "-----------\n",join(":",@ksnp),"\n";
		if (" @ksnp "=~/ 0 1 2 3 /) { 
			splice @ksnp,0,4;
			$trim=4;
			}
		elsif (" @ksnp "=~/ 0 1 3 /) { 
			splice @ksnp,0,3;
			$trim=4;
			}
		elsif (" @ksnp "=~/ 0 2 3 /) { 
			splice @ksnp,0,3;
			$trim=4;
			}
		elsif (" @ksnp "=~/ 0 1 2 /) { 
			splice @ksnp,0,3;
			$trim=3;
			}
		elsif (" @ksnp "=~/ 0 1 /) { 
			splice @ksnp,0,2;
			$trim=2;
			}
		elsif (" @ksnp "=~/ 1 2 /) { 
			splice @ksnp,0,2;
			$trim=3;
			}
		if (" @ksnp "=~/ 32 33 34 35 /) { 
			splice @ksnp,-4;
			$trim=-4;
			}
		elsif (" @ksnp "=~/ 32 34 35 /) { 
			splice @ksnp,-3;
			$trim=-4;
			}
		elsif (" @ksnp "=~/ 32 33 35 /) { 
			splice @ksnp,-3;
			$trim=-4;
			}
		elsif (" @ksnp "=~/ 33 34 35 /) { 
			splice @ksnp,-3;
			$trim=-3;
			}
		elsif (" @ksnp "=~/ 34 35 /) { 
			splice @ksnp,-2;
			$trim=-2;
			}
		elsif (" @ksnp "=~/ 33 34 /) { 
			splice @ksnp,-2;
			$trim=-3;
			}
#if ($trim!=0) { print STDERR join(":",@ksnp),"\ntrim:$trim\n";}
	}
	if (!@ksnp) {
#print STDERR "no snps left\n";
		if ($monokeep==1) { 
			$goodsite{$r}=1;
			$tagprint++;
			next;
			}
		else { next; }
	}
	@{$snps{$r}}=@ksnp;
	if ($trim>0) { 
		foreach (my $sp=0; my $spos=${$snps{$r}}[$sp];$sp++) { ${$snps{$r}}[$sp]=$spos-$trim+1;}
#print STDERR join(":",@{$snps{$r}}),"\n";
	}
	else { 
		foreach ($sp=0; $spos=${$snps{$r}}[$sp];$sp++) { ${$snps{$r}}[$sp]=$spos+1;}
#if ($trim!=0) { print STDERR join(":",@{$snps{$r}}),"\n";}
	}
	$trimtag{$r}=$trim;

	$goodsite{$r}=1;
	$tagprint++;
	
	if ($#ksnp>0) { 
#print STDERR "phasing\n";
		$phase++;
		$isphase=$phase;
	}
	else { $isphase=0;}
	foreach my $sn (@ksnp){
		my @freqs=();
		next if ($sn=~/HASH/);
		my @bases=split ("",$snp{$sn});
#print "\tSNP at $sn (@bases):\n";
		my %gt1={};
		my %gt2={};
		my $nsam=0;
		my $skew=0;
		my $het=0;
		foreach $s (@samples) {
			next if (!${$gta{$s}}[0]);
			$nsam++;
			my @a1=split("",${$gta{$s}}[0]);
			my @a2;
			if (${$gta{$s}}[1]){
				@a2=split("",${$gta{$s}}[1]);
				$skew=$skew + (${$seen2{$s}}{${$gta{$s}}[1]}/(${$seen2{$s}}{${$gta{$s}}[1]}+${$seen2{$s}}{${$gta{$s}}[0]}));
				$het++;
#print "\tsample $s:\n\t",@a1,"\n\t",@a2,"\n";
			}
			else { 
				@a2=@a1;
#print "\tsample $s: \n\t",@a1," homo\n";
			
			}
			for ($b=0; $ba=$bases[$b];$b++) {
				if ($a1[$sn] eq $ba ) { 
#print "\t\t\tb $b : allele1\t$sn\t$ba\n";
					$gt1{$s}=$b; 
					$freqs[$b]++;
				}
				if ($a2[$sn] eq $ba ) { 
#print "\t\t\tb $b : allele2\t$sn\t$ba\n";
					$gt2{$s}=$b; 
					$freqs[$b]++;
				}
			}
		}
		if ($het) {$skew=$skew/$het;}
		else {$skew=1;}
		if ($skew>1) { $skew=1/$skew;}
		$skew=sprintf("%.2f",$skew);
		$skew=$skew*100;
		my $rdep=$sumcount;
		for ($b=0; $ba=$bases[$b];$b++) {
			$freqs[$b]=$freqs[$b]/(2*$nsam);
			$freqs[$b]=sprintf("%.3f",$freqs[$b]);
		}
#print "\t\tcoord $coord skew $skew; rdep $rdep; freqs ",join(':',@freqs),"\n";
		$tagg=$r;
		print {TAB} "$chrom\t",$sn+1,"\t.\t",shift @bases,"\t",join(",",@bases),"\t$skew\t.\tSB=$sbr;AB=$skew;NS=$nsam;TP=",$sn+1,";DP=$rdep;AF=",join(":",@freqs),"\tGT:GQ:PS:AD:DP";
		foreach $s (@samples){
#print "\t\t\tSAMPLE $s, genotype $gt1{$s}|$gt2{$s} isphase? $isphase ($phase)\n";			
			if (${$gta{$s}}[0]) {
				my $seca;
				my $td;
				my $gq;
				if (${$seen2{$s}}{${$gta{$s}}[1]}) { 
					$seca=${$seen2{$s}}{${$gta{$s}}[1]};
					$td=${$seen2{$s}}{${$gta{$s}}[0]}+$seca;
					if ($td>40) { $gq=40; } else { $gq=$td;}
#print "${$gta{$s}}[0]:${$seen2{$s}}{${$gta{$s}}[0]},${$gta{$s}}[1]:$seca:TD=$td\n";
					$gq=sprintf("%.0f",20+(-10)*(log(2*(0.5**$gq))/log(10)));
				}
				else { 
					$seca =".";
					$td=${$seen2{$s}}{${$gta{$s}}[0]};
					if ($td<$mindp) { goto NOCALL;}
					if ($td>40) { $gq=40; } else { $gq=$td;}
					$gq=sprintf("%.0f",(-10)*(log(2*(0.5**$gq))/log(10)));
				}
				if ($isphase) {
#					if (!${$seen2{$s}}{${$gta{$s}}[1]}) { $seca=${$seen2{$s}}{${$gta{$s}}[0]};}
					print {TAB} "\t$gt1{$s}|$gt2{$s}:$gq:$phase:",${$seen2{$s}}{${$gta{$s}}[0]},",",$seca,":",$td;
				}
				else {
					print {TAB} "\t$gt1{$s}/$gt2{$s}:$gq:.:",${$seen2{$s}}{${$gta{$s}}[0]},",",$seca,":",$td;
				}
			}
			else {
				NOCALL:
				print {TAB} "\t./.:.:.:.,.:.";
			}												
		}
		print {TAB} "\n";
	}
#	undef $goodmapQ{$r};	
}
close TAB;

$now=localtime;
warn "Haplo-VCF writing start: $now\n";

# writing VCF-haplo format
$tabname=$gmap."_Vhap_count".$cut."_ab".$minQ."_sb".$strb."_clip".$clip.".vcf";
open TAB, ">$tabname";
print {TAB} "##fileformat=VCFv4.1
##INFO=<ID=SB,Number=1,Type=Integer,Description=\"strand bias (ratio x100)\">
##INFO=<ID=AB,Number=1,Type=Integer,Description=\"minor allele bias (ratio x100)\">
##INFO=<ID=Trim,Number=1,Type=Integer,Description=\"the first starting or last remaining base after trimming\">
##INFO=<ID=TP,Number=1,Type=Integer,Description=\"SNP positions in the tag\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"combined read depth\">
##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled genotype quality\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele depth\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
foreach $sa (@samples) {
	print {TAB} "\t$sa";
}
print {TAB} "\n";
my $chrom;
my $pos;
foreach $r (sort {$a <=> $b} @goodrefs){
#print "ref $r\n";
	next if ($r=~/HASH/);
	next if (!$goodsite{$r});
	$chrom=$r;
	my $maj="";
	my $maxcount=0;
	my $sumcount=0;
	foreach $ss (@{$goodmapQ{$r}}) {
		if ($seen{$ss}>$maxcount){
			$maxcount=$seen{$ss};
			$maj=$ss;
		}
		$sumcount+=$seen{$ss};
	}
	my @ordmap=();	
	push @ordmap, $maj;
	
	my $sbr=0;
	my %trimmed={};
	foreach $ss (@{$goodmapQ{$r}}) {
		$sbr+=$sb{$ss};
		my $sss;
#		my $majs;
		if ($trimtag{$r}>=0 ) {
			$sss=substr($ss,$trimtag{$r});
#			$majs=substr($maj,$trimtag{$r});
		}
		else {
			$sss=substr($ss,0,36+$trimtag{$r});
#			$majs=substr($maj,0,36+$trimtag{$r});
		}
		$trimmed{$ss}=$sss;
		next if ($ss eq $maj);
		push @ordmap, $ss;
	}
	next if ($#ordmap==0 && $monokeep==0);
	$sbr=sprintf("%.0f",$sbr/($#{$goodmapQ{$r}}+1));

	my %gt={};
	my %gta={};
	my %acounts={};
	my %anum={};
	my @ordmapt;
	foreach $s (@samples) {
		for ($a=0;$ssa=$ordmap[$a];$a++) {
			if (" @ordmapt "!~/ $trimmed{$ssa} /){
				push @ordmapt,$trimmed{$ssa};
				$anum{$trimmed{$ssa}}=$#ordmapt;
			}
			if(${$seen2{$s}}{$ssa}) {				
#				push @{$gt{$s}}, $a;
				$acounts{$s}{$trimmed{$ssa}}+=${$seen2{$s}}{$ssa};
				if (" @{$gta{$s}} "!~/$trimmed{$ssa}/) {
					push @{$gta{$s}},$trimmed{$ssa};
				}
			}
		}
	}
#print join("\n",@ordmap),"\n";
	my $skew=0;
	my $nsam=0;
	my $het=0;
	my %freq={};
	foreach $s (@samples){
#print "\t\t\tsample $s, genotype ",join('/',@{$gt{$s}}), " @{$gta{$s}}\n";
		if ($#{$gta{$s}}>1) { die "more than two alleles at $r $ss\n";}
		elsif ($#{$gta{$s}}>0) {
			$skew=$skew + ($acounts{$s}{${$gta{$s}}[1]}/($acounts{$s}{${$gta{$s}}[1]}+$acounts{$s}{${$gta{$s}}[0]}));
			$het++;
			$freq{${$gta{$s}}[1]}++;
			$freq{${$gta{$s}}[0]}++;
			$nsam++;
		}
		elsif (${$gta{$s}}[0]) { 
			$nsam++; 
			$freq{${$gta{$s}}[0]}+=2;
		}
	}
	if ($het) {$skew=$skew/$het;}
	else {$skew=1;}
	if ($skew>1) { $skew=1/$skew;}
	$skew=sprintf("%.0f",100*$skew);
	my $rdep=$sumcount;
	my @freqs;
	foreach $all (keys %freq ){
		next if ($all=~/HASH/);
		$freq{$all}=$freq{$all}/(2*$nsam);
		$freq{$all}=sprintf("%.3f",$freq{$all});
		push @freqs,$freq{$all} ;
	}
	if ($#ordmapt==0) { 
		push @ordmapt, "."; 
		$trimtag{$r}=0;
		push @{$snps{$r}},100;
	}
	$maj=shift @ordmapt;
#print "\t\t\tnsam=$nsam skew=$skew het=$het freqs=@freqs\n";
#	if ($trimtag{$r}>=0) {
#		$maj=substr($maj,$trimtag{$r});
#		for (my $ssa=0;$ssa<=$#ordmap;$ssa++) { $ordmap[$ssa]=substr($ordmap[$ssa],$trimtag{$r});}
#	}
#	else { 
#		$maj=substr($maj,0,36+$trimtag{$r});
#		for ($ssa=0;$ssa<=$#ordmap;$ssa++) { $ordmap[$ssa]=substr($ordmap[$ssa],0,36+$trimtag{$r});}
#	}
	print {TAB} "$chrom\t1\t.\t$maj\t",join(",",@ordmapt),"\t$skew\t.\tSB=$sbr;AB=$skew;Trim=$trimtag{$r};NS=$nsam;DP=$rdep;TP=",join(",",@{$snps{$r}}),";AF=",join(":",@freqs),"\tGT:GQ:AD:DP";
	foreach $s (@samples){
		my $seca;		
		my $td;
		my $gq;
		if ($acounts{$s}{${$gta{$s}}[1]}) { 
			$seca=$acounts{$s}{${$gta{$s}}[1]};
			$td=$acounts{$s}{${$gta{$s}}[0]}+$seca;
#print "${$gta{$s}}[0]:${$seen2{$s}}{${$gta{$s}}[0]},${$gta{$s}}[1]:$seca:TD=$td\n";
			if ($td>40) { $gq=40; } else { $gq=$td;}
			$gq=sprintf("%.0f",20+(-10)*(log(2*(0.5**$gq))/log(10)));
		}
		else { 
			$seca =".";
			$td=$acounts{$s}{${$gta{$s}}[0]};
			if ($td<$mindp) { goto NOCALL2;}
			if ($td>40) { $gq=40; } else { $gq=$td;}
			$gq=sprintf("%.0f",(-10)*(log(2*(0.5**$gq))/log(10)));
		}
		if ($#{$gta{$s}}>0) {
#print "\t\t\tSAMPLE $s, genotype ",join('/',@{$gt{$s}}), " @{$gta{$s}}\n";			
			print {TAB} "\t",$anum{${$gta{$s}}[0]},"/",$anum{${$gta{$s}}[1]},":",$gq,":",$acounts{$s}{${$gta{$s}}[0]},",",$seca,":",$td;
		}
		elsif (${$gta{$s}}[0]) {
			print {TAB} "\t",$anum{${$gta{$s}}[0]},"/",$anum{${$gta{$s}}[0]},":",$gq,":",$acounts{$s}{${$gta{$s}}[0]},",",$seca,":",$td;
		}				
		else {
			NOCALL2:
			print {TAB} "\t./.:.:.,.:.";
		}												
	}
	print {TAB} "\n";
	undef $goodmapQ{$r};
}
close TAB;
$now=localtime;
warn "all done: $now\n";




