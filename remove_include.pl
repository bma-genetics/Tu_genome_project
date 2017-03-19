#!/usr/bin/perl
#Author: huilong du
#Note: 去除包含关系
use warnings;
use strict;
use Bio::SeqIO;
use Bio::Perl;


my $infile =  shift;          # blasr selected result   03
my $outfile = shift;          # remove included ctg     
open IN,"<$infile" or die $!;
open OUT,">$outfile" or die $!;
my %existed_ctg=();

while(<IN>){
    chomp;
    my $line=$_;
    my @content=split "\t",$line;
    next if($content[0] eq "include");
    print OUT "$line\n";
}
close(IN);
close(OUT);
