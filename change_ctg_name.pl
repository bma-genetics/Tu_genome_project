#!/usr/bin/perl
#Author: huilong du
#Note: changing the ctg name for non-repeated
use warnings;
use strict;
my $infile=shift;
my $outfile=shift;

open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;

my @content=split "/",$infile;
$content[-1]=~/(\S+)-bac.fasta/;
my $sign_bac=$1;

my $count=1;
while(<IN>){
     chomp;
     my $line=$_;
     if($line=~/^>/){
         print OUT ">$sign_bac-$count\n";
         $count++;
     }
     else{
         print OUT "$line\n";
     }
}
close(IN);
close(OUT);
