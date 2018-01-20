#!/usr/bin/perl

my $usage= "

Prints out list of commands for launcher_creator.py 
to bowtie2-map the trimmed 36b 2bRAD reads

Arguments:
1: glob to fastq files
2: reference to map to (basename of bowtie2 index file)
3: optional, the position of name-deriving string in the file name
	if separated by underscores, 
	such as: input file Sample_RNA_2DVH_L002_R1.cat.fastq
	specifying arg2 as \'3\' would create output file with a name \'2DVH.fastq'	
			
";

if (!$ARGV[0]) { die $usage;}
my $glob=$ARGV[0];
if (!$ARGV[1]) { die $usage;}
my $ref=$ARGV[1];


opendir THIS, ".";
my @fqs=grep /$glob/,readdir THIS;
my $outname="";

foreach $fqf (@fqs) {
	if ($ARGV[2]) {
		my @parts=split('_',$fqf);
		$outname=$parts[$ARGV[1]-1].".bt2.sam";
	}
	else { $outname=$fqf.".bt2.sam";}
	print "bowtie2 --no-unal --score-min L,16,1 --local -L 16 -x $ref -U $fqf -S $outname\n";
}
