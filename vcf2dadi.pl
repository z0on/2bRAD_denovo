#!/usr/bin/env perl
use strict;
use warnings;
#use Bio::SeqIO;

my $usage="

Usage: perl vcf2dadi.pl <vcf file> <list file>

Genome file is in fasta format;
List file gives population designations, like this:

sample1    population1
sample2    population1
sample3    population2

";

if (!$ARGV[1]) { die $usage;}
my ($VCF,$list)=@ARGV;

my %list;
open(IN,"< $list")||die"$!";
while (<IN>) {
    chomp;
    next if(/^\#/);
    my @a=split(/\s+/);
    $list{$a[0]}=$a[1];
}
close IN;

my %vcf;
my %pop;
my %record;
my %crec;
my @chroms;
open(IN,"< $VCF");
while (<IN>) {
    chomp;
    next if(/^##/);
    if(/^#/){
        my @a=split(/\s+/);
        for(my $i=9;$i<@a;$i++){
            next if(!exists $list{$a[$i]});
#            die"$a[$i] does not exists in $list\n" if(!exists $list{$a[$i]});
            $record{$i}=$list{$a[$i]};
        }
        next;
    }
    my @a=split(/\s+/);
    my ($chr,$pos,$ref,$alt)=($a[0],$a[1],$a[3],$a[4]);
    if (!$crec{$chr}) {
	    push @chroms, $chr;
	    $crec{$chr}=1;
#print "@chroms\n-----\n";
	}
    next if($alt=~/,/);

    $vcf{$chr}{$pos}{ref}=$ref;
    $vcf{$chr}{$pos}{alt}=$alt;
#   print "--------------\n$ref\t$alt\n";
    foreach my $i(keys %record){
        my $indv=$record{$i};
        $pop{$indv}=1;
#       print "$indv\t";
        my $geno=$a[$i];
#       print "$geno\t";
        $geno=~/^(.)[\/\|](.):/;
        my $tempa=$1;
        my $tempb=$2;
#       print "$tempa $tempb\n";
        my ($a1,$a2)=(0,0);
        if($tempa eq "." || $tempb eq "."){
            $a1=0;
            $a2=0;
        }
        elsif ($tempa+$tempb == 1) {
            $a1=1;
            $a2=1;
        }
        elsif ($tempa+$tempb == 2) {
            $a1=0;
            $a2=2;
        }
        elsif ($tempa+$tempb == 0) {
            $a1=2;
            $a2=0;
        }
        $vcf{$chr}{$pos}{$indv}{a1}+=$a1;
        $vcf{$chr}{$pos}{$indv}{a2}+=$a2;
    }
    #last;
}
close IN;

$VCF=~s/\..+//;
my $outname=$VCF."_dadi.data";
open(O,"> $outname");
my $title="NAME\tOUT\tAllele1";
foreach my $pop(sort keys %pop){
    $title.="\t$pop";
}
$title.="\tAllele2";
foreach my $pop(sort keys %pop){
    $title.="\t$pop";
}
$title.="\tGene\tPostion\n";
print O "$title";
foreach my $id (@chroms){
#print "chrom $id...\n";
    foreach my $pos (sort {$a<=>$b} keys %{$vcf{$id}}){
        my $ref="A".$vcf{$id}{$pos}{ref}."A";
        my $line="$ref\t$ref\t$vcf{$id}{$pos}{ref}";
        foreach my $pop(sort keys %pop){
            my $num=$vcf{$id}{$pos}{$pop}{a1};
            $line.="\t$num";
        }
        $line.="\t$vcf{$id}{$pos}{alt}";
        foreach my $pop(sort keys %pop){
            my $num=$vcf{$id}{$pos}{$pop}{a2};
            $line.="\t$num";
        }
        $line.="\t$id\t$pos\n";
        print O "$line";
    }
}
close O;

