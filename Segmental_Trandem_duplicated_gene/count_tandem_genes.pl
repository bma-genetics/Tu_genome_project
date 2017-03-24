#!/usr/bin/perl -w
use strict;
##count tandem duplicated genes in Turartu
my @tempArr;
my $input;

my $output1 = "count_tandem_genes.out";
open(OUTFILE1, ">".$output1) ||
    die "Unable to open $output1 for writing: $!\n";

##read gff file to save information on location of each gene from Tu.
for(my $i=1;$i<=7;$i++){
    $input = "/mnt/hgfs/D/Wheat_Agenome_DNA/Scan_dup_block/TuvsTu_duplication/tandem_genes/Tu".$i."_tandem_genes.out";
    open(INFILE, $input) ||  
	die "Unable to open $input for reading: $!\n";
    
    my %tandem_counter = ();

    while(<INFILE>){                   
	chomp;
	if($_ =~ /^The/){next;}
	@tempArr = split(/\s+/,$_);
	$tandem_counter{$tempArr[0]}= "";
    }
    close INFILE;
    my $counter = 0;
    for my $key (keys %tandem_counter){
	$counter++;
    }
    print OUTFILE1 "Tu".$i."\t".$counter."\n";
    close INFILE;
}
close OUTFILE1;
