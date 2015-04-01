#!/usr/bin/perl

my $usage="

vhap2migrate.pl :

converts \"Vhap\" vcf file (produced by haplocall_denovo) to MIGRATE-n sequence input
(make sure you keep monomorphic loci when running haplocall_denovo: mono=keep)

Arguments: 

           arg1 : input vhap vcf file
arg2 (optional) : comma-delimited list of files listing samples from specific  
                  populations

Output: 
MIGRATE-n formatted genotypes (printed to STDOUT)

Example: 
vhap2migrate.pl denovo.vhap.filt.recode.vcf keppel.pop,orpheus.pop,maggi.pop 

";

$infile=shift or die $usage;
open VCF, $infile or die "cannot open input file $infile\n";

my @popfiles=();
if ($ARGV[0]) { @popfiles=split(",",$ARGV[0]);}

my %sampop={};
my %popsam={};
my @pops=();

if ($#popfiles>0){
	foreach my $pop (@popfiles) {
		open POP, $pop or die "cannot open pop file $pop\n";
		$pop=~s/\..+//;
		while (<POP>){
			chomp;
			push @{$sampop{$pop}},$_;
			$popsam{$_}=$pop;
		}
		push @pops, $pop;
	}
}

my %gtypes={};
my @loci;
my %llen={};
my %inds={};

while (<VCF>) {
	next if ($_=~/^##/);
	if ($_=~/^#CHROM/) {
		chomp;
		@samples=split("\t",$_);
		splice @samples, 0, 9;
		next;
	}
	chomp;
	my @dline=split ("\t",$_);
	my $locus=$dline[0]."_".$dline[1];
	push @loci, $locus;
	$length{$locus}=length($ref);
	my $ref=$dline[3];
	$llen{$locus}=length($ref);
	my @alleles=split(",",$dline[4]);
	unshift @alleles, $ref;
	for(my $i=0;$samples[$i];$i++){
		next if (!$popsam{$samples[$i]});
		my @gt=split(":",$dline[9+$i]);
		@gta=split("/",$gt[0]);
		if ($gta[0] ne ".") {
			push @{$gtypes{$locus}{$popsam{$samples[$i]}}}, ($alleles[$gta[0]],$alleles[$gta[1]]);
			push @{$inds{$locus}{$popsam{$samples[$i]}}},("$samples[$i]_1","$samples[$i]_2");
		}
	}
}

print "N ",$#pops+1," ",$#loci+1,"\n";
my @ls;
foreach my $locus (@loci){
	push @ls, $llen{$locus};
}
print join(" ",@ls),"\n";

foreach my $pop (@pops) {
	my @ls=();
	foreach $locus (@loci) {
		push @ls, $#{$inds{$locus}{$pop}}+1;
	}
	print join(" ",@ls)," $pop\n";
	foreach $locus (@loci) {
		for ($i=0;$i<=$#{$inds{$locus}{$pop}};$i++){
			print ${$inds{$locus}{$pop}}[$i],"\t",${$gtypes{$locus}{$pop}}[$i],"\n";
		}
	}
}