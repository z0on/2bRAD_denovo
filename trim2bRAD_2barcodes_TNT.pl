#!/usr/bin/perl
$usage="

2bRAD-TNT reads trimmer (May 2014)

Trims 2b-RAD fastq files to leave only those with adaptor on the far end and restriction site;

Deduplicating: removes all but one read sharing the same 3-base degeneracy within the TNT
adaptor, the same first 34 bases of the read, and the same secondary barcode

arguments:
input=[fastq filename]
site=[perl-style pattern for restrictase site] 
		Default is BcgI: \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\'
bc2=[perl-style barcode pattern] in-line barcode that immediately follows the RAD fragment, 
		default \'[ATGC]{4}\'
adaptor=[far-end adaptor sequence] default ATCG (as in TNT adaptor)
clip=[integer] number of bases to clip off the ends of the reads (corresponding to 
       ligated overhangs), default 0
minBCcount=[integer] minimum count per in-line barcode to output a separate file. 
				     Default 10000.
sampleID=[integer] the position of name-deriving string in the file name
					if separated by underscores, such as: 
					for input file Sample_RNA_2DVH_L002_R1.cat.fastq
					specifying arg2 as \'3\' would create output 
					file with a name \'2DVH.trim'	
";

my $input="";
if ("@ARGV"=~/input=(\S+)/){ $input=$1; }
else { die $usage; }
open INP, $1 or die "cannot open input $1\n";
my $site='.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}';
if ("@ARGV"=~/site=(\S+)/){ $site=$1;}
my $bcod="[ATGC]{4}";
if ("@ARGV"=~/bc2=(\S+)/){ $bcod=$1;}
my $degen="AG[ATGC]A[ATGC]A[ATGC]";
my $adap="ATCG";
if ("@ARGV"=~/adaptor=(\S+)/){ $adap=$1;}
my $clip=0;
if ("@ARGV"=~/clip=(\d+)/){ $clip=$1;}
my $sampleid=100;
if("@ARGV"=~/sampleID=(\d+)/){ $sampleid=$1;}
my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";
my $minBCcount=10000;
if ("@ARGV"=~/minBCcount=(\d+)/){$minBCcount=$1;}

my $ll=3;
my %data={};
my %dups={};
my $counter=0;
my %dupnum={};
while (<INP>) {
	chomp;
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		$counter++;
#print "$seq:";
		if ($seq=~/^($site)($bcod)($degen)$adap/) {
#print "$1:$2:$3\n";
			$chk=$1.$2.$3;
			if ($dups{$chk}) {
				$dupnum{$2}++;
#print "dup found: $chk\n";
			}
			else {
				$dups{$chk}=1;
				my $rd=substr($1,$clip,length($1)-$clip*2);
				$qua=substr($qua,$clip,length($1)-$clip*2);
				$dline="$name bcd=$2\n$rd\n+\n$qua\n";
				push @{$data{$2}}, $dline ;
			}
#print "$dline\n";
		}
		else {
#print "\n\nNOT FINDING $site:$bcod:$adap\n\n";
		}
#print "-----------\n$counter\n$_\n";
		$seq="";
		$ll=0;
		$qua="";
		@sites=();
		$name=$name2;
	}
	elsif ($ll==0){
		chomp;
		$seq=$_;
		$ll=1;
	}
	elsif ($ll==2) { 
		chomp;
		$qua=substr($_,3,length($_)); 
		$ll=3;
	}
	else { $ll=2;}
}
$name2=$1; 
if ($seq=~/^($site)($bcod)($degen)$adap/) {
#	my $chk=substr($2,0,34);
	$chk=$1.$2.$3;
	if ($dups{$chk}) {
		$dupnum{$1}++;
#print "dup found: $chk\n";
	}
	else {
		$dups{$chk}=1;
		my $rd=substr($1,$clip,length($1)-$clip*2);
		$qua=substr($qua,$clip,length($1)-$clip*2);
		$dline="$name bcd=$2\n$rd\n+\n$qua\n";
		push @{$data{$2}}, $dline ;
	}
}

my $outname;
if ($sampleid<100) {
	my @parts=split('_',$input);
	$outname=$parts[$sampleid-1];
}
else { 
	$input=~s/\..+//;
	$outname=$input;
}

my $totgood=0;
my $totdups=0;
foreach $bc2 (sort keys %data) {
	next if ($bc2=~/HASH/);
	next if ($#{$data{$bc2}}<$minBCcount);
	print $outname,"_","$bc2: goods: ",$#{$data{$bc2}}+1," ; dups: $dupnum{$bc2}\n";
	$totgood=$totgood+$#{$data{$bc2}}+1;
	$totdups=$totdups+$dupnum{$bc2};
	my $outname1=$outname."_".$bc2.".tr0";
	open OUT, ">$outname1" or die "cannot create $outname1\n";
	foreach $d (@{$data{$bc2}}) { print {OUT} $d; }
	close OUT;
}

warn "\n$outname: total goods : $totgood ; total dups : $totdups\n";