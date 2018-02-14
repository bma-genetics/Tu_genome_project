#Usage:calculate average expression within sliding windows along chromsomes
use strict;
use warnings;

my $infile;
my $outfile;
my @tempArr;
my %exp_dic = ();
my %long_gene = ();
my %temp_dic = ();
my $tissue;
my @tempArr1;

##read the expression data
$infile = "/mnt/hgfs/D/Wheat_Agenome_mRNA/A_trans_reads/FPKM_DE_3tissues.out";
open INFILE, "$infile" or die "Could not open $infile: $! \n";

while(my $line = <INFILE>){
    chomp($line);
    if($line =~ /^Transcript/){next;}
    @tempArr = split(/\t/,$line);
    $tissue = "leaf";
    $temp_dic{$tissue}{$tempArr[0]}=(log2($tempArr[2]+1)+log2($tempArr[3]+1)+log2($tempArr[4]+1)+log2($tempArr[5]+1))/4;
    $tissue = "spike";
    $temp_dic{$tissue}{$tempArr[0]}=(log2($tempArr[6]+1)+log2($tempArr[7]+1)+log2($tempArr[8]+1))/3;
    $tissue = "root";
    $temp_dic{$tissue}{$tempArr[0]}=(log2($tempArr[9]+1)+log2($tempArr[10]+1))/2;
    @tempArr1 = split(/\./,$tempArr[0]);
    $long_gene{$tempArr1[0]}{$tempArr[1]} = $tempArr[0];
}
close(INFILE);

for my $key (sort keys %long_gene){
    for my $subkey (sort {$b<=>$a} keys %{$long_gene{$key}}){
	for $tissue (sort keys %temp_dic){
	    $exp_dic{$tissue}{$key} = $temp_dic{$tissue}{$long_gene{$key}{$subkey}};
	}
	last;
    }
}

%temp_dic = ();
##read location info of each gene
$infile = "/mnt/hgfs/D/Wheat_Agenome_DNA/Agenome.v3.0/WheatTu_Chr.IGDBV1_Final_gene/WheatTu_Chr.IGDBV1_Final.gff";
open INFILE, "$infile" or die "Could not open $infile: $! \n";

my %loc_dic = ();
while(my $line = <INFILE>){
    chomp($line);
    @tempArr = split(/\s+/,$line);
    if($tempArr[2] =~ /gene/){
	my $temp = $tempArr[3] + int(($tempArr[4] - $tempArr[3] + 1)/2);
	my @tempArr1 = split(/\./,$tempArr[8]);
	$loc_dic{substr($tempArr[0],2)}{$temp} = substr($tempArr1[0],3);
    }
}
close(INFILE);

#for my $key (sort {$a<=>$b} keys %loc_dic){
#    for my $subkey (sort {$a<=>$b} keys %{$loc_dic{$key}}){
#	print $key."\t".$subkey."\t".$loc_dic{$key}{$subkey}."\n";
#    }
#}exit;


##calculate average expression level in each sliding windows
my $start = 0;
my $end = 1;
my $win_size = 5000000;
my %total_count = ();
my %exp_count = ();
my $chr_len;

for(my $i = 1;$i <=7;$i++){
    if($i eq 1){
	$chr_len = 570000000;
    }
    elsif($i eq 2){
	$chr_len = 741000000;
    }
    elsif($i eq 3){
	$chr_len = 734000000;
    }
    elsif($i eq 4){
	$chr_len = 605000000;
    }
    elsif($i eq 5){
	$chr_len = 650000000;
    }
    elsif($i eq 6){
	$chr_len = 572000000;
    }
    else{
	$chr_len = 713000000;
    }
    $outfile = "chr".$i."_exp_total.out";
    open (OUTFILE, ">$outfile") or die "Cannot open $outfile for writing \n";
    print OUTFILE "Win_start\tWin_end\tgeneNumb\tExp_leaf\troot\tspike\n";
    
    %total_count = ();
    %exp_count = ();
    
    for($start = 1;$start < $chr_len;$start+=$win_size){
	$end = $start+$win_size;
	for my $key (sort {$a<=>$b} keys %{$loc_dic{$i}}){
	    if($key < $start){
		delete $loc_dic{$i}{$key};next;
	    }
	    if($key > $end){last;}
	    for(my $j= 0;$j <3;$j++){
		if($j==0){$tissue = "leaf";}
		elsif($j==1){$tissue = "spike";}
		else{$tissue = "root";}
		if(!exists $total_count{$start}{$tissue}){
		    $total_count{$start}{$tissue} = 0;
		    $exp_count{$start}{$tissue} = 0;
		}
		if(exists $exp_dic{$tissue}{$loc_dic{$i}{$key}}){
		    $total_count{$start}{$tissue}++;
		    $exp_count{$start}{$tissue}+= $exp_dic{$tissue}{$loc_dic{$i}{$key}};
		}
	    }
	}
    }
    for my $key (sort {$a<=>$b}  keys %exp_count){
	print OUTFILE $key."\t".($key + $win_size)."\t";
	for my $subkey (sort keys %{$exp_count{$key}}){
	    #print $subkey."\n";
	    my $temp; 
	    if($subkey eq "leaf"){
		print OUTFILE $total_count{$key}{$subkey}."\t";
		#print $total_count{$key}{$subkey}."\t";
	    }
	    if($total_count{$key}{$subkey} >0){
		#$temp = $exp_count{$key}{$subkey}/$total_count{$key}{$subkey};
		$temp = $exp_count{$key}{$subkey};
	    }	    
	    else{$temp = 0;}
	    print OUTFILE $temp."\t";
	}
	print OUTFILE "\n";
    }
    close OUTFILE;
}

print "The running of sliding_win_exp_density.pl is done\n";

sub log2 {
my $n = shift;
return log($n)/log(2);
}
