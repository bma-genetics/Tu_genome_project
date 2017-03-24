#!/usr/bin/perl -w
use strict;
##pick tandem duplicated genes within Turartu
my @tempArr;
my $input;

my %segmental_counter = ();
##read duplicated blocks
for(my $i=1;$i<=7;$i++){
    for(my $j=$i;$j<=7;$j++){
	$input = "/mnt/hgfs/D/Wheat_Agenome_DNA/Scan_dup_block/TuvsTu_duplication/segmental_genes/Tu".$i."_Tu".$j.".blocks";
	open(INFILE, $input) ||  
	    die "Unable to open $input for reading: $!\n";
	
	while(<INFILE>){                   
	    chomp;
	    if($_ =~ /TuG1812/){
		@tempArr = split(/\s+/,$_);
		$segmental_counter{$i}{$tempArr[0]}= "";
		$segmental_counter{$j}{$tempArr[2]}= "";
	    }
	}
	close INFILE;
	$input = "/mnt/hgfs/D/Wheat_Agenome_DNA/Scan_dup_block/TuvsTu_duplication/segmental_genes/Tu".$i."_Tu".$j."_reverse.blocks";
	open(INFILE, $input) ||  
	    die "Unable to open $input for reading: $!\n";
	
	while(<INFILE>){                   
	    chomp;
	    if($_ =~ /TuG1812/){
		@tempArr = split(/\s+/,$_);
		$segmental_counter{$i}{$tempArr[0]}= "";
		$segmental_counter{$j}{$tempArr[2]}= "";
	    }
	}
	close INFILE;
    }
}

#my $output1 = "count_segmental_genes.out";
#open(OUTFILE1, ">".$output1) ||
#    die "Unable to open $output1 for writing: $!\n";

##read segmental genes
for(my $i = 1;$i <= 7;$i++){
    my $counter = 0;
    for my $key (keys %{$segmental_counter{$i}}){
	print $i."\t".$key."\n";
	$counter++;
    }
    print OUTFILE1 "Tu".$i."\t".$counter."\n";
    close INFILE;
}
close OUTFILE1;
