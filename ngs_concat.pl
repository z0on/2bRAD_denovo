#!/usr/bin/perl

$usage= "

concatenates files by matching pattern in their names

arg1: common pattern for files
arg2: perl-like pattern in the filename to recognize, 
	  use brackets to specify the unique part, as in
	\"FilenameTextImmediatelyBeforeSampleID(.+)FilenameTextImmediatelyAfterSampleID\"

Example (to concatenate files names like Sample_Pop1_L1.fastq, Sample_Pop1_L2.fastq): 
ngs_concat.pl 'Sample' 'Sample_(.+)_L'

";

my $ff = shift or die $usage;
my $patt=shift or die $usage;
#print "pattern $patt\n";
opendir THIS, ".";
my @files=grep /$ff/,readdir THIS;
print "files:\n",@files, "\n";
my @ids=();
foreach $file (@files){
	if ($file=~/$patt/) {
		$ii=$1;
#		unless (grep {$_ eq $ii} @ids){ 
#			my $name=$file;
#			my $ccat=$patt;
#			$ccat=~s/\(.+\)/$ii/;
#			$ccat.="*";
			$name=$ii.".fq";
print "$file > $name\n";

			`cat $file >> $name`;
#			push @ids, $ii;
#		}
	}
	else { print "$patt not found in $file\n";}
}

