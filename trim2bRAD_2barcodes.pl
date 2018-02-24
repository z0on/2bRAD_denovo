#!/usr/bin/perl
$usage="

Trims 2b-RAD fastq files to leave only those with 
adaptor on the far end and restriction site

            fastq=[file name] : fastq filename
    site=[perl-style pattern] : pattern defining the site, default
                                \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\' for BcgI
barcode2=[perl-style pattern] : in-read barcode that immediately follows the RAD fragment, 
                                 default \'[ATGC]{4}\'
adaptor=[sequence or pattern] : adaptor sequence, default AGATCGGAAG
               trim=[integer] : number of bases to trim off the ends of the reads 
                                default 0
         minBCcount=[integer] : minimum count per in-line barcode to output 
                                a separate file. Default 10000.
           sampleID=[integer] : the position of name-deriving string in the file name 
                                if separated by underscores, such as: 
			            		for input file Sample_RNA_2DVH_L002_R1.cat.fastq
			            		specifying arg2 as \'3\' would create output 
		           		    	file with a name \'2DVH.trim'	
		           		    	
";
my $infile="";
if ("@ARGV"=~/fastq=(\S+)/) { 
	open INP, $1 or die $usage; 
	$infile=$1;
}
else { die $usage;}
my $site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}";
if ("@ARGV"=~/site=(\S+)/){$site=$1;}
my $bcod="[ATGC]{4}";
if ("@ARGV"=~/barcode2=(\S+)/){$bcod=$1;}
my $adap="AGATCGGAAG";
if ("@ARGV"=~/adaptor=(\S+)/){$adap=$1;}
my $clip=0;
if ("@ARGV"=~/trim=(\d+)/){$clip=$1;}
my $sampleid=100;
if ("@ARGV"=~/sampleID=(\d+)/){$sampleid=$1;}
my $minBCcount=10000;
if ("@ARGV"=~/minBCcount=(\d+)/){$minBCcount=$1;}
my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";

my $ll=3;
my %data={};
my $counter=0;
while (<INP>) {
	chomp;
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		$counter++;
#print "$seq:";
		if ($seq=~/^($site)($bcod)$adap/) {
#print "$1:$2\n";
			my $rd=substr($1,$len,length($1)-$clip*2);
			$qua=substr($qua,$clip,length($1)-$clip*2);
			$dline="$name bcd=$2\n$rd\n+\n$qua\n";
			push @{$data{$2}}, $dline ;
#print "$dline\n";
		}
		else {
#print "\t\tNOT FINDING $site:$bcod:$adap\n";
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
		$qua=substr($_,length("$,length($_)); 
		$ll=3;
	}
	else { $ll=2;}
}
$name2=$1; 
if ($seq=~/^($site)($bcod)$adap/) {
	my $rd=substr($1,$clip,length($1)-$clip*2);
	$qua=substr($qua,$clip,length($1)-$clip*2);
	$dline="$name bcd=$2\n$rd\n+\n$qua\n";
	push @{$data{$2}}, $dline ;
}

my $outname;
if ($sampleid<100) {
	my @parts=split('_',$infile);
	$outname=$parts[$sampleid-1];
}
else { 
	$infile=~s/\..+//;
	$outname=$infile;
}

#print "\n\nOUT: $outname\n";

foreach $bc2 (sort keys %data) {
	next if ($bc2=~/HASH/);
#print "BC: $bc2  COUNT: ",$#{$data{$bc2}}," MIN: $minBCcount\n"; 
	next if ($#{$data{$bc2}}<$minBCcount);
	my $outname1=$outname."_".$bc2.".tr0";
	open OUT, ">$outname1" or die "cannot create $outname1\n";
	foreach $d (@{$data{$bc2}}) { print {OUT} $d; }
	close OUT;
}
