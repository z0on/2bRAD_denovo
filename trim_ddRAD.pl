#!/usr/bin/env perl
$usage="

Trims ddRAD R1 fastq files for analysis using Agnostic Protocol

            fastq=[file name] : fastq filename
    site=[perl-style pattern] : pattern defining the site, default \'CATGC\' for SphI
 barcode=[perl-style pattern] : in-read barcode that immediately preceding the restriction site, 
                                 default \'[ATGCN]{5}\'
     excise=[integer-integer] : bases to excise from the middle of the read (1-based). Default 0-0.
             length=[integer] : final length of read after trimming; 
                                 default 100
         minBCcount=[integer] : minimum count per in-line barcode to output 
                                a separate file. Default 100000.
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
my $site="CATGC";
if ("@ARGV"=~/site=(\S+)/){$site=$1;}
my $bcod="[ATGCN]{5}";
if ("@ARGV"=~/barcode=(\S+)/){$bcod=$1;}
my $ex1=0;
my $ex2=0;
if ("@ARGV"=~/excise=(\d+)-(\d+)/){
	$ex1=$1;
	$ex2=$2;
}
#print "excising: $ex1-$ex2\n\n";
my $sampleid=100;
if ("@ARGV"=~/sampleID=(\d+)/){$sampleid=$1;}
my $minBCcount=10000;
if ("@ARGV"=~/minBCcount=(\d+)/){$minBCcount=$1;}
my $len=100;
if ("@ARGV"=~/length=(\d+)/){$len=$1;}

my $trim=0;
my $name="";
my $name2="";
my $seq="";
my $qua="";
my $dummy="";
my $ll=3;
my %data={};
my $counter=0;
while (<INP>) {
	chomp;
	if ($ll==3 && $_=~/^(\@.+)$/) {
		$name2=$1; 
		$counter++;
#print "$seq\n";
		if ($ex1>0 && $seq) { 
			my $seq1=substr $seq, 0,$ex1-1; 
			my $seq2=substr $seq, $ex2,length($seq)-1;
			$seq=$seq1.$seq2; 
#print "$seq\n";
		}
		if ($seq=~/^($bcod)($site)(.{$len})/) {
#print "$1:$2:$3\n";
			my $rd=$3;
			$qua=substr($qua,length("$1$2"),$len);
			$dline="$name bcd=$1\n$rd\n+\n$qua\n";
			push @{$data{$1}}, $dline ;
#print "\n$dline\n";
		}
		else {
#print "\t\tNOT FINDING $bcod:$site \n";
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
		$qua=$_; 
		if ($ex1>0) { 
			my $seq1=substr $qua, 0,$ex1-1; 
			my $seq2=substr $qua, $ex2,length($seq)-1;
			$qua=$seq1.$seq2; 
		}
		$ll=3;
	}
	else { $ll=2;}
}
$name2=$1; 
$name2=$1; 
$counter++;
#print "$seq\n";
if ($ex1>0 && $seq) { 
	my $seq1=substr $seq, 0,$ex1-1; 
	my $seq2=substr $seq, $ex2,length($seq)-1;
	$seq=$seq1.$seq2; 
#print "$seq\n";
}
if ($seq=~/^($bcod)($site)(.{$len})/) {
#print "$1:$2:$3\n";
	my $rd=$3; 
	$qua=substr($qua,length("$1$2")-1,$len);
	$dline="$name bcd=$1\n$rd\n+\n$qua\n";
	push @{$data{$1}}, $dline ;
#print "\n$dline\n";
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