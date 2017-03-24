#!/usr/bin/perl -w
use strict;

#distance <= one gene
#block with at least two genes

my @tempArr;
my @tempArr1;
my @tempArr2;
my $input;

##read gff file to save information on location of each gene from Tu.
$input = "/mnt/hgfs/D/Wheat_Agenome_DNA/Agenome.v3.0/WheatTuGene.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %Tu_id = ();
my %Tu_sort = ();
my $sort = 0;
while(<INFILE>){                   
   chomp;
   @tempArr = split(/\s+/,$_);
   if($tempArr[2] eq "gene"){
       @tempArr1 = split(/[=;]/,$tempArr[8]);
       $sort++;
       $Tu_sort{$tempArr1[1]} = $sort;
   }
}
close INFILE;

for(my $i = 1;$i<=7;$i++){
    my $output1 = "Tu".$i."_tandem_genes.out";
    open(OUTFILE1, ">".$output1) ||
	die "Unable to open $output1 for writing: $!\n";
  
    ##scan tandem duplicated blocks on chromosome  i
    my %homo_group = ();
    my $group_num = 0;
    my $marker = 0;
    $input = "/mnt/hgfs/D/Wheat_Agenome_DNA/Scan_dup_block/TuvsTu_duplication/Tu".$i."vsTu".$i."_pair.out";
    open(INFILE, $input) ||  
	die "Unable to open $input for reading: $!\n";
   
    while(<INFILE>){    
	chomp;
	@tempArr = split(/\s+/,$_);
	$marker = 0;
	for (my $j=1;$j<= $group_num;$j++){
	    if(exists $homo_group{$j}{$tempArr[0]} || exists $homo_group{$j}{$tempArr[3]}){
		$homo_group{$j}{$tempArr[0]} = "";
		$homo_group{$j}{$tempArr[3]} = "";
		$marker = 1;
		last;
	    }
	}
	if($marker == 0){
	    $group_num ++;	
	    $homo_group{$group_num}{$tempArr[0]} = "";
	    $homo_group{$group_num}{$tempArr[3]} = "";
	}
    }
    close INFILE;
    #for (my $j=1;$j<= $group_num;$j++){
    #	print "genes in group ".$j.":\n";
    #	for my $id (keys %{$homo_group{$j}}){
    #	    print $id."_".$Tu_sort{$id}."\t";
    #	}
    # print "\n";
    #}exit;

    my $block_num = 0;
    for (my $j=1;$j<= $group_num;$j++){
	my $ext = 0;
	my %temp_block = ();
	my %sorted_group = ();
	for my $id (keys %{$homo_group{$j}}){
	    $sorted_group{$Tu_sort{$id}} = $id;
	}
	my @gene_dic = ();
	for my $key (sort {$a<=>$b} keys %sorted_group){
	    push @gene_dic, $sorted_group{$key};
	}
	#for(my $j=0;$j<= $#gene_dic;$j++){
	#    print $gene_dic[$j]."_".$Tu_sort{$gene_dic[$j]}."\t";
	#}exit;
	for(my $j=0;$j< $#gene_dic;$j++){
	    if($Tu_sort{$gene_dic[$j+1]} - $Tu_sort{$gene_dic[$j]} <=2){
		if($ext == 0){
		    $temp_block{$Tu_sort{$gene_dic[$j]}} = $gene_dic[$j];
		    $ext++;
		}
		$temp_block{$Tu_sort{$gene_dic[$j+1]}} = $gene_dic[$j+1];
		$ext++;
	    }
	    elsif($ext >=2){
		$block_num++;
		print OUTFILE1 "The tandem duplicate block number:".$block_num." on Tu".$i."\n";
		for my $key (sort  {$a<=>$b} keys %temp_block){
		    print OUTFILE1 $temp_block{$key}."\t".$key."\n";
		}
		$ext = 0;
		%temp_block = ();
	    }
	}
	if($ext >=2){
	    $block_num++;
	    print OUTFILE1 "The tandem duplicate block number:".$block_num." on Tu".$i."\n";
	    for my $key (sort  {$a<=>$b} keys %temp_block){
		print OUTFILE1 $temp_block{$key}."\t".$key."\n";
	    }
	    $ext = 0;
	    %temp_block = ();
	}
    }
    close OUTFILE1;
}  
 


