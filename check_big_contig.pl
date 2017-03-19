#!/usr/bin/perl
#Author: huilong du
#Note: 检查大的，没有被连起来的contig 
use warnings;
use strict;
my $infile=shift;               #list all non-connected contigs file
open IN,"<$infile" or die $!;
my $outfile=shift;              #list all bac >5k non-connected contigs
open OUT,">$outfile" or die $!;
my $big_contig_size=shift;

my %existed_bac=();
while(<IN>){
    chomp;
    my $line=$_;
    my @field=split "/",$line;
    my $bac_name=$field[-2];
    open BAC,"<$line" or die $!;
    my $sign="";
    while(<BAC>){
        chomp;
        my $contig=$_;
        if($contig=~/^>(\S+)/){
             $sign=$1;
        }
        else{
             my $len=length($contig);
             if($len>$big_contig_size && $sign=~/contig_size/){
                    print OUT "$bac_name\t$sign\t$len\n";
                    $existed_bac{$bac_name}++;
             }
        }
    }
    close(BAC);
}
foreach my $key (keys %existed_bac){
    print "$key\t$existed_bac{$key}\n";
}
close(OUT);
                   


