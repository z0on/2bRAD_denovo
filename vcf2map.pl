#!/usr/bin/env perl

my $usage="

vcf2map.pl : replaces scaffold coordinates in the vcf file with their approximate locations 
in the linkage map. Scaffolds not pinned by the map are discarded.
The input VCF file is assumed to be sorted by scaffold.

Arguments:

vcf=[filename]
map=[filename] : table where third column is linkage group, fourth column is map position,
                 fifth column is the scaffold ID, and sixth column is scaffold position.
       cmmb=3  : recombination rate (centimorgans per megabase). Set to 0 to estimate 
                 from data on a per-scaffold basis (be sure to check with 
                 MapAnchoring_tests.R that the results are reasonable)
mindist=10000  : Minimum base distance to calculate recombination rate (cM/Mb)

Outputs three files mapAnchored_[input vcf filename] :
.vcf: recoded sorted vcf
.rec: inferred ecombination rate across genome, cM/Mb
.tab: table of new to old coordinates:
        chromosome	chrom.coordinate	scaffold	scaff.coordinate	
.rec: table of inferred recombination rates (cM/Mb - tends to be super noisy though): 
        linkageGroup	scaffoldSegment	meanCoordinate	recombRate

Example: 
vcf2map.pl vcf=gatkFinal2digitifera.vcf map=oldMap.tab   

Mikhail Matz, matz\@utexas.edu
May 2015

";

# to calculate means:
use List::Util qw(sum);
sub mean {
    return sum(@_)/@_;
}


my $vcf="";
my $mapfile="";
my $cmmb=3;
#my $rescale=0;
my $mindist=2500;
if ("@ARGV"=~/vcf=(\S+)/) { $vcf=$1; } else { die $usage;}
if ("@ARGV"=~/map=(\S+)/) { $mapfile=$1; } else { die $usage;}
if ("@ARGV"=~/cmmb=(\S+)/) { $cmmb=$1; }
if ("@ARGV"=~/mindist=(\S+)/) { $mindist=$1; }
#if ("@ARGV"=~/rescale/) { $rescale=1; }

open MAP, $mapfile or die "cannot open map table $map\n";
my %map={};
my %lg={};
my %coord={};
my %lgcount={};
my @lgs=();
while (<MAP>){
	chop;
	my @line=split("\t",$_);
	if (" @{$map{$line[4]}} "!~/$line[3]/) {
		push @{$lg{$line[4]}},$line[2];
#warn "lgs for $line[4]: @{$lg{$line[4]}}\n";
		push @lgs, $line[2] unless (" @lgs "=~/ $line[2] /);
		$lgcount{$line[4]}{$line[2]}++;
		push @{$map{$line[4]}}, $line[3] ;
		push @{$coord{$line[4]}}, $line[5];
	}
}
close MAP;

my @scaffolds=keys %lg;

# extracting consistently mapped scaffold segments:
# for multi-marker scaffolds, these are regions with at least two consecutive markers 
# mapping to the same linkage group
my %onelg={};
my $onemap={};
my $onecoord={};
my %splits={};
my %lg2scaff={};
my $singletons=0;
foreach my $s (@scaffolds){
	next if ($s=~/HASH/);
	if (@{$lg{$s}}==1) {
		my $sgl=$s.${$lg{$s}}[0];
		$onelg{$sgl}=${$lg{$s}}[0];
		push @{$lg2scaff{$l}},$sgl;
		push @{$onemap{$sgl}},${$map{$s}}[0];
		push @{$onecoord{$sgl}},${$coord{$s}}[0];
		push @{$splits{$s}},$sgl;
		$singletons++;
	}
	else {	
		my @goodlgs=();
		for(my $i=0; my $l=${$lg{$s}}[$i]; $i++){
			if ($l == ${$lg{$s}}[$i-1]) { 
				push @goodlgs,$l unless (" @goodlgs "=~/ $l /);
			}
		}
		if (@goodlgs==0) {
#warn "\t\tambiguous mapping for $s: @{$lg{$s}}\n";
				next;
		}
		foreach my $gl (@goodlgs) {
			my $sgl=$s.$gl;
			push @{$splits{$s}},$sgl;
			for(my $i=0; my $l=${$lg{$s}}[$i]; $i++){
				if ($l eq $gl) {				
					$onelg{$sgl}=$l;
					push @{$lg2scaff{$l}},$sgl;
					push @{$onemap{$sgl}},${$map{$s}}[$i];
					push @{$onecoord{$sgl}},${$coord{$s}}[$i];
				}
			}
#warn "$sgl: @{$lg{$s}} => $onelg{$sgl}: @{$onecoord{$sgl}}\n";
		}
	}
}
@scaffolds=keys %onelg;
warn "\n",$#scaffolds+1," scaffolds\n",$singletons," of them are singletons\n";


# calculating recombination rates and scaffold orientations
# in multi-marker scaffolds, require monotonously increasing marker base coordinates 
my %direction={};
my %meanmap={};
my @allrecs=();
foreach my $s (@scaffolds){
	next if ($s=~/HASH/);
#warn "\n---------------------\n$s:  @{$onemap{$s}}  :  @{$onecoord{$s}}\n";
	if (@{$onemap{$s}}==1) {
		$direction{$s}=1;
		$meanmap{$s}=${$onemap{$s}}[0];
	}
	else {
		my $mid=$#{$onemap{$s}}/2;
		my $beg=0;
		my $end=0;
		my $all=0;
		my $i=0;
		my $monotone=1;
		my $prevchange=0;
		my @ratios;
		for ($i=0;$i<@{$onemap{$s}};$i++){ 
			if ($i>0 and $monotone) {
				my $change=${$onecoord{$s}}[$i]-${$onecoord{$s}}[$i-1];		
				my $changem=${$onemap{$s}}[$i]-${$onemap{$s}}[$i-1];		
				if (abs($change)>$mindist){
					if ($change*$prevchange<0 ) { 
						$monotone=0;
#warn "\tSHIFT $i change:$change;prev:$prevchange mono:$monotone\n";
					}
					else { 
#warn "\t\t$i change:$change;prev:$prevchange prod:",$change*$prevchange," mono:$monotone\n";
						$prevchange=$change;
						push @ratios,abs($changem/$change);
					}
				}
			}
			$meanmap{$s}+=${$onemap{$s}}[$i];
			if ($i<$mid) {$beg+=${$onecoord{$s}}[$i];}
			elsif($i>$mid) {$end+=${$onecoord{$s}}[$i];}
		}
		$meanmap{$s}=sprintf("%.4f",$meanmap{$s}/@{$onemap{$s}});
		if ($beg<=$end){$direction{$s}=1;}
		else {$direction{$s}=-1;}
		if ($prevchange and $monotone) {
			$rec{$s}=mean(@ratios)*1e+6;
			push @allrecs,@ratios;
#warn ">>>monotonous: @ratios MEAN:",mean(@ratios),"\n";			 
		}
	}
#warn "\t\tmean:$meanmap{$s} direction:$direction{$s}\n";
}
my $gr=mean(@allrecs)*1e+6;
@allrecs=sort{$a<=>$b} @allrecs;
my $grmed=$allrecs[sprintf("%.0f",@allrecs/2)]*1e+6;
warn "\nmean cM/Mb:\t$gr\nmedian cM/Mb:\t$grmed\n(",scalar(@allrecs)," comparisons)\n\n";
if ($cmmb!=0) { $grmed=$cmmb;}
warn "using cM/Mb = $grmed\n";


my $outvcf="mapAnchored_".$vcf;
open OUTV, ">$outvcf" or die "cannot create output vcf file $outvcf\n";
open VCF, $vcf or die "cannot open vcf file $vcf\n";

# re
my %data={};
my %tab={};
my $chr="";
my $add=0;
my $seen={};
while (<VCF>) {
	if ($_=~/^#/) { print {OUTV} $_ and next;}
	(my $s,my $co,my @rest)=split(/\t/,$_);	
	next unless @{$splits{$s}}>0;
	my $nearest=1e+9;
	my $spot=0;
	for (my $i=0; my $sgl=${$splits{$s}}[$i];$i++){
		my $c=mean(@{$onecoord{$sgl}});
		if (abs($co-$c)<$nearest) {
			$split=$sgl;
			$nearest=abs($co-$c);
		}
	}
	my $scale=$grmed;
#	if ($rec{$split} and $rescale){
#warn "using native rec.rate for $split: $rec{$split}\n";
#		$scale=$rec{$split};
#	}
#if (!$direction{$s}) { warn "no dir: $s  $meanmap{$split}  @{$onecoord{$split}}\n";}
	my $loc=sprintf("%.0f",1e+6*$meanmap{$split}/$scale+$direction{$split}*($co-sprintf("%.0f",mean(@{$onecoord{$split}}))));
#warn "$s: co:$co; meanmap:$meanmap{$split}; meanCoord:",mean(@{$onecoord{$split}}),"; loc:$loc\n"; 
	my $ind=1e+10*$onelg{$split}+$loc;
	if ($seen{$ind}) { 
warn "$split   $co   repeated index:$ind   mm:$meanmap{$split}    @{$onecoord{$split}}   $loc\n$tab{$ind}\n";
	$loc+=1;
	$ind+=1;
	}
	if ($seen{$ind}) { 
warn "\t$split   $co   repeated index:$ind   mm:$meanmap{$split}    @{$onecoord{$split}}   $loc\n\t$tab{$ind}\n";
	$loc+=1;
	$ind+=1;
	}
	$seen{$ind}++;
	$data{$ind}=join("\t",@rest);
	$data{$ind}="chr$onelg{$split}\t$loc\t".$data{$ind};
	$tab{$ind}="chr$onelg{$split}\t$loc\t$s\t$co\n";
#warn "$s:$co:  $onelg{$split}    @{$onecoord{$split}}     $loc\n";
}

my @indx=sort {$a <=> $b } keys %data;
my $count=0;
my $chr="";
my $add=0;
my $ch="";
my $l=0;
my @rest=();
foreach my $i (@indx) {  
	next if ($i=~/HASH/);
	($ch,$l,@rest)=split(/\t/,$data{$i});
	if ($ch ne $chr){
		$chr=$ch;
		$add=0-$l;
#warn "new chrom: $ch   add: $add\n";
	}
	$l+=$add;
	$data{$i}=$ch."\t".$l."\t".join("\t",@rest);
	(my $cht, my $lt,my @restt)=split("\t",$tab{$i});
	$tab{$i}=$ch."\t".$l."\t".join("\t",@restt);
	$count++;
}
warn "\n$count variants anchored\n\n";

my $outtab=$outvcf;
$outtab=~s/\.vcf/\.tab/;
open OUTT, ">$outtab" or die "cannot create output tab file $outtab\n";

foreach my $i (@indx) {  
	next if ($i=~/HASH/);
	print {OUTV} $data{$i};
	print {OUTT} $tab{$i};
}
close OUTV;
	
my $outtab=$outvcf;
$outtab=~s/\.vcf/\.rec/;
open OUTT, ">$outtab" or die "cannot create output tab file $outtab\n";
foreach my $lg (@lgs){
	my %seen={};
	my @scafs =sort { $meanmap{$a} <=> $meanmap{$b} } @{$lg2scaff{$lg}};
#warn "lg$lg   @scafs\n";
	foreach my $s (@scafs){
		next if ($s=~/HASH/ or $seen{$s} or !$rec{$s});
		print {OUTT} "$lg\t$s\t",sprintf("%.0f",1e+6*$meanmap{$s}/$cmmb),"\t$rec{$s}\n";
		$seen{$s}=1;
	}
}
	


