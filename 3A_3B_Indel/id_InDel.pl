#!/usr/bin/perl -w
use strict;

##identify insertion and deletion of genes in Ta3B and Tu3
##Using bd/Os as reference to decide whether a gene insertion/deletion have ever occurred
##Only mind genes have no syntenic orthologs between 3A and 3B
##If an absent gene exists at least two species, the gene was defined to DELETION
##If 1) exist genes have no homologs in all other species and 2) have reference genes at both bounder and 3) connected genes are not more than 10. 

my @tempArr;
my @tempArr1;
my $input;
my $output;


#set key as reference start position; 0 set as no indel ; 1 is indel
my @in_dic = ();
my @del_dic = ();

##mark the closest line with both Tu/Ta genes or another is defined as deletion.
##contains two elements the first is sort of the line and the second is status 
##status: 0 is not bd/Os homologs; 1 is have a bd or os homologs; 2 is both bd and os homologs
my $before_pos = 0;
my $before_status = 0;
##pos counter start from 0
my $pos_counter = -1;

my @block_dic = ();
my $block_id;

for(my $sp = 0;$sp < 2;$sp++){
    @in_dic = ();
    @del_dic = ();
    @block_dic = ();
    $before_pos = 0;
    $before_status = 0;
    $pos_counter = -1;
    $block_id = "";
    ##read Tu_ref_collinearity.out to identify Ta deletion and Tu insertion
    if($sp == 0){
	$input = "Ta_ref_collinearity.out";
	$output = "Ta_ref_collinearity_indel.out";
    }
    else{
	$input = "Tu_ref_collinearity.out";
	$output = "Tu_ref_collinearity_indel.out";
    }
    open(INFILE, $input) ||  
	die "Unable to open $input for reading: $!\n";
    open(OUTFILE, ">".$output) ||
	die "Unable to open $output for writing: $!\n";

    if($sp == 0){
	print OUTFILE "TaID\tTastart\tTaend\tInsertion\tTuID\tTustart\tTuend\tDeletion\tBdID\tBdstart\tBdend\tOsID\tOsstart\tOsend\tSbID\tSbstart\tSbend\n";
    }
    else{
	print OUTFILE "TuID\tTustart\tTuend\tInsertion\tTaID\tTastart\tTaend\tDeletion\tBdID\tBdstart\tBdend\tOsID\tOsstart\tOsend\tSbID\tSbstart\tSbend\n";
    }
    while(<INFILE>){    
	chomp;
	if($_ =~ /^tuID/ || $_ =~ /^taID/){next;}
	if($_ =~ /Alignment/){
	    if($pos_counter > 0){
		print OUTFILE $block_id."\n";
		for(my $j=$before_pos;$j<=$pos_counter;$j++){
		    push(@in_dic,0);
		}
		$pos_counter = -1;
		for (my $line =0;$line <=$#block_dic;$line++){
		    $pos_counter++;
		    @tempArr = split(/\s+/,$block_dic[$line]);
		    print OUTFILE $tempArr[0]."\t".$tempArr[1]."\t".$tempArr[2]."\t".$in_dic[$pos_counter];
		    print OUTFILE "\t".$tempArr[3]."\t".$tempArr[4]."\t".$tempArr[5]."\t".$del_dic[$pos_counter];
		    print OUTFILE "\t". $tempArr[6]."\t".$tempArr[7]."\t".$tempArr[8]."\t". $tempArr[9]."\t".$tempArr[10]."\t".$tempArr[11]."\t". $tempArr[12]."\t".$tempArr[13]."\t".$tempArr[14]."\n";
		}
	    } 
	    @in_dic = ();
	    @del_dic = ();
	    @block_dic = ();
	    $before_pos = 0;
	    $before_status = 0;
	    $pos_counter = -1;
	    $block_id = $_;
	    next;
	}
	$pos_counter ++;
	push(@block_dic,$_);
	@tempArr = split(/\s+/,$_);
	if($tempArr[3] eq "--"){
	    ##2 or 3 of Bd/Os/Sb are there
	    if($tempArr[6] ne "--" && $tempArr[9] ne "--" || $tempArr[6] ne "--" && $tempArr[12] ne "--" || $tempArr[9] ne "--" && $tempArr[12] ne "--"|| $tempArr[6] ne "--" && $tempArr[9] ne "--" && $tempArr[12] ne "--" ){
		push(@del_dic,1);
		##search insertion from before to previous line
		search_insertion();
		$before_pos = $pos_counter;
		$before_status = 2;
	    } 
	    ##only one in there
	    elsif($tempArr[6] ne "--" || $tempArr[9] ne "--"|| $tempArr[12] ne "--"){
		push(@del_dic,0);
		##search insertion from before to previous line
		#search_insertion();
		for(my $i=$before_pos;$i<$pos_counter;$i++){
		    push(@in_dic,0);
		}
		$before_pos = $pos_counter;
		$before_status = 1;
	    }
	    else{
		push(@del_dic,0);
	    }
	}
	else{      
	    ##search insertion from before to previous line
	    push(@del_dic,0);
	    ##2 or 3 of Bd/Os/Sb are there
	    if($tempArr[6] ne "--" && $tempArr[9] ne "--"|| $tempArr[6] ne "--" && $tempArr[12] ne "--" || $tempArr[9] ne "--" && $tempArr[12] ne "--"|| $tempArr[6] ne "--" &&$tempArr[9] ne "--" && $tempArr[12] ne "--"){
		search_insertion();
		$before_status = 2;
	    } 
	    ##only one there
	    elsif($tempArr[6] ne "--" || $tempArr[9] ne "--"|| $tempArr[12] ne "--"){
		#search_insertion();
		for(my $i=$before_pos;$i<$pos_counter;$i++){
		    push(@in_dic,0);
		}
		$before_status = 1;
	    }
	    else{
		for(my $i=$before_pos;$i<$pos_counter;$i++){
		    push(@in_dic,0);
		}
		$before_status = 0;
	    }
	    $before_pos = $pos_counter;
	}
    }
    close INFILE;

    for(my $j=$before_pos;$j<=$pos_counter;$j++){
	push(@in_dic,0);
    }

    #for(my $j=0;$j<=$#in_dic;$j++){
    #    print $j."\t".$in_dic[$j]."\n";
    #}
    #exit;
    ##pos counter start from 0
    $pos_counter = -1;
    print OUTFILE $block_id."\n";
    for(my $line =0;$line <= $#block_dic;$line++){
	$pos_counter++;
	@tempArr = split(/\s+/,$block_dic[$line]);
	print OUTFILE $tempArr[0]."\t".$tempArr[1]."\t".$tempArr[2]."\t".$in_dic[$pos_counter];
	print OUTFILE "\t".$tempArr[3]."\t".$tempArr[4]."\t".$tempArr[5]."\t".$del_dic[$pos_counter];
	print OUTFILE "\t". $tempArr[6]."\t".$tempArr[7]."\t".$tempArr[8]."\t".$tempArr[9]."\t".$tempArr[10]."\t".$tempArr[11]."\t".$tempArr[12]."\t".$tempArr[13]."\t".$tempArr[14]."\n";
    }
    close OUTFILE;
}


sub search_insertion{
    ##set before line
    if($pos_counter > 0){
	push(@in_dic,0);
    }
    ##set lines after before line
    #print $before_status."\t".$before_pos."\t".$pos_counter."\n";
    if($before_status > 1){
	for (my $i= $before_pos+1;$i< $pos_counter;$i++){
	    push(@in_dic,1);
	}
    }
    else{
	for (my $i= $before_pos+1;$i< $pos_counter;$i++){
	    push(@in_dic,0);
	}
    }
}

print "end of id_InDel.pl\n";
