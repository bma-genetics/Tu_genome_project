#!/usr/bin/perl -w
use strict;

##count insertion and deletion for collinear Ta and Tu genes.

my @tempArr;
my @tempArr1;
my $input;
my $output;

#count genes with consensus insertion and deletion estimation
my %in_dic = ();
my %del_dic = ();
my %delete_in = ();
my %delete_del = ();

for(my $sp = 0;$sp < 2;$sp++){
    %in_dic = ();
    %del_dic = ();
    %delete_in = ();
    %delete_del = ();
    my $in_count = 0;
    my $del_count = 0;
    ##read Tu_ref_collinearity.out to identify Ta deletion and Tu insertion
    if($sp == 0){
	$input = "Ta_ref_collinearity_indel.out";
    }
    else{
	$input = "Tu_ref_collinearity_indel.out";
    }
    open(INFILE, $input) ||  
	die "Unable to open $input for reading: $!\n";
    open(OUTFILE, ">".$output) ||
	die "Unable to open $output for writing: $!\n";

    while(<INFILE>){    
	chomp;
	if($_ =~ /^TuID/ || $_ =~ /^TaID/ || $_ =~ /Alignment/){next;}
	@tempArr = split(/\s+/,$_);
	if(exists $in_dic{$tempArr[0]}){
	    if($in_dic{$tempArr[0]} ne $tempArr[3]){
		delete($in_dic{$tempArr[0]});
		$delete_in{$tempArr[0]}= "";
	    }
	}
	elsif(exists $delete_in{$tempArr[0]}){
	    ;
	}
	else{
	    $in_dic{$tempArr[0]}=$tempArr[3];
	}
	if(exists $del_dic{$tempArr[0]}){
	    if($del_dic{$tempArr[0]} ne $tempArr[7]){
		delete($del_dic{$tempArr[0]});
		$delete_del{$tempArr[0]}= "";
	    }
	}
	elsif(exists $delete_del{$tempArr[0]}){
	    ;
	}
	else{
	    $del_dic{$tempArr[0]}=$tempArr[7];
	}
    }
    close INFILE;
    close OUTFILE;
    for my $key (keys %in_dic){
	if($in_dic{$key} > 0){
	    $in_count++;
	}
    }
    for my $key (keys %del_dic){
	if($del_dic{$key} > 0){
	    $del_count++;
	}
    }
    print $in_count."\n";
    print $del_count."\n";
}



print "end of count_InDel.pl\n";
