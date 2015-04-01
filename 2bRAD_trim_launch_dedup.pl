#!/usr/bin/perl

my $usage= "

2bRAD_trim_launch_dedup.pl :

Prints out list of calls to trim2bRAD_2barcodes_dedup.pl, one for
each fastq file, to to trim and deduplicate raw 2bRAD reads.

See help for trim2bRAD_2barcodes_dedup.pl for more details.

prints to STDOUT

Arguments:
arg1, required: glob to fastq files
site=[pattern] perl-style pattern of the restriction site to recognise, in single quotes.
		default is BcgI: \'.{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}\'
		if you need AlfI, here is how to define it: \'.{12}GCA.{6}TGC.{12}\'
adaptor=[dna sequence] adaptor sequence to look for on the far end of the read. 
		Default AGATC
sampleID=[integer] the position of name-deriving string in the file name
					if separated by underscores, such as: 
					for input file Sample_RNA_2DVH_L002_R1.cat.fastq
					specifying arg2 as \'3\' would create output 
					file with a name \'2DVH.trim'	
barcode2=[integer] length of the in-line barcode immediately following 
			   the restriction fragment. Default 0. 
			   
Example:
2bRAD_trim_launch_dedup.pl fastq sampleID=3 > trims
	
";

my $glob=shift or die $usage;
my $site=".{12}CGA.{6}TGC.{12}|.{12}GCA.{6}TCG.{12}";
if("@ARGV"=~/site=(\S+)/){ $site=$1;}
my $sampleid=100;
if("@ARGV"=~/sampleID=(\d)/){ $sampleid=$1;}
my $adaptor="AGATC";
if("@ARGV"=~/adaptor=(\S+)/){ $adaptor=$1;}


opendir THIS, ".";
my @fqs=grep /$glob/,readdir THIS;
my $outname="";
my @outnames;

foreach $fqf (@fqs) {
	if ($sampleid<100) {
		my @parts=split('_',$fqf);
		$outname=$parts[$sampleid-1].".tr0";
	}
	else { $outname=$fqf.".tr0";}
	print "trim2bRAD_2barcodes_dedup.pl input=$fqf site=\"$site\" adaptor=\"$adaptor\" sampleID=$sampleid\n";
}
