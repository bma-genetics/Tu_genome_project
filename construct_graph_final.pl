#!/usr/bin/perl
use strict;
no warnings;
use Graph;
use Graph::Directed;
use Graph::Undirected;
use Graph::Easy;
use Bio::SeqIO;
use Bio::Perl;
use Data::Dumper;

my $DEBUG=0;
my $inputfile=shift;                #blasr_result output
my $ContigNameListFile=shift;       #all contig name       contig_name.list
my $MergedContigFile = shift;       #output graph png 
#my $out_cluster_file= shift;        #output ctg clusters
my @ColumnNameArray = qw ( RefPos RefStrand Identity Score RefName QryName QryPos QryStart QryEnd QryLength RefStart RefEnd RefLength );

my ( $ContigGraph , $EasyGraph ) = &initializeGraph( $ContigNameListFile );
my ( $Finding_Best_Graph , $Output_Graph) = &initializeGraph( $ContigNameListFile );

my %NodeHash = %{ ( &GenerateRecordHash( $inputfile ) )[0] };
my %EdgeHash = %{ ( &GenerateRecordHash( $inputfile  ))[1] };
my %NodePairsHash = %{ ( &GenerateNodePairHash( \%EdgeHash ))[0] };
#&drawEasyGraph ( $EasyGraph , $MergedContigFile , '01' );


############ 整合所有节点和边的信息，合并边，并且选择最好的pacbio作为其连接，并且计算两个NODE连接边的权重 #######
my %ctg_ctg_count=%{ (&count_ctg_pair( \%NodePairsHash))[0]};
my %ctg_ctg_overlap_len=%{ (&count_ctg_pair(\%NodePairsHash))[1] };
my %ctg_ctg_ori=%{ (&count_ctg_pair(\%NodePairsHash))[2]};
my %ctg_ctg_chain=%{ (&count_ctg_pair(\%NodePairsHash))[3]};
my %ctg_ctg_line=%{ (&count_ctg_pair(\%NodePairsHash))[4]};


=pod

############ 过滤图中的边，如果一个节点的度>=3，则选择不同方向的两条边 ##########################################
%ctg_ctg_count = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[0])};
%ctg_ctg_overlap_len=%{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[1])};
%ctg_ctg_ori = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[2])};
%ctg_ctg_chain = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[3])};
%ctg_ctg_line = %{ (&Filtering_Edge (\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line))[4])};
########### 结束 ################################################################################################
=cut

############ 根据节点和边，将每对contig分别进行延伸，每次选择最好的边  ##########################################
my @ctg_clusters = &Finding_best_pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line, \%NodePairsHash );
#print "@ctg_clusters\n";

&draw_ctg_clusters(\@ctg_clusters);


if($DEBUG==1){
    print Dumper(%NodePairsHash);
    print Dumper(%ctg_ctg_count);
    print Dumper(%ctg_ctg_overlap_len);
    print Dumper(%ctg_ctg_ori);
    print Dumper(%ctg_ctg_chain);
    print Dumper(%ctg_ctg_line);
}






############## 统计支持ctg对的pacbio数，以及选择overlap最长的pacbio作为其连接两个ctg的信息 ########################
sub count_ctg_pair{
    my $CtgPairsHash = shift;
    my %NodePair=%{ $CtgPairsHash };
    my %ctg_ctg_count=();
    my %ctg_ctg_overlap_len=();
    my %ctg_ctg_ori=();
    my %ctg_ctg_chain=();
    my %ctg_ctg_line=();
    open PAIR,">ctg_pairs.txt" or die $!;
    open ORI,">ctg_ctg_ori.txt" or die $!;
    foreach my $key1 (keys %NodePair){
        my %CtgPairPacbio= %{ $NodePair{$key1} };
        my @content=split "-",$key1;
        my $left_ctg=$content[0];
        my $right_ctg=$content[1];
        my $best_pacbio="";
        foreach my $pacbio (keys %CtgPairPacbio){
             $ctg_ctg_count{$left_ctg}{$right_ctg}++;
             my $overlap_len=$CtgPairPacbio{$pacbio}{'First'}{'RefEnd'}-$CtgPairPacbio{$pacbio}{'First'}{'RefStart'}+$CtgPairPacbio{$pacbio}{'Next'}{'RefEnd'}-$CtgPairPacbio{$pacbio}{'Next'}{'RefStart'};
             if(!exists $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}){
                     $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}=0;
             }
             if($ctg_ctg_overlap_len{$left_ctg}{$right_ctg}<$overlap_len){
                     $ctg_ctg_overlap_len{$left_ctg}{$right_ctg}=$overlap_len;
                     $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}=$CtgPairPacbio{$pacbio}{'First'}{'RefPos'};
                     $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}=$CtgPairPacbio{$pacbio}{'Next'}{'RefPos'};
                     $ctg_ctg_ori{$right_ctg}{$left_ctg}{1}=$CtgPairPacbio{$pacbio}{'Next'}{'RefPos'};
                     $ctg_ctg_ori{$right_ctg}{$left_ctg}{2}=$CtgPairPacbio{$pacbio}{'First'}{'RefPos'};
                     
                     $ctg_ctg_chain{$left_ctg}{$right_ctg}{1}=$CtgPairPacbio{$pacbio}{'First'}{'RefStrand'};
                     $ctg_ctg_chain{$left_ctg}{$right_ctg}{2}=$CtgPairPacbio{$pacbio}{'Next'}{'RefStrand'};
                     $ctg_ctg_chain{$right_ctg}{$left_ctg}{1}=$CtgPairPacbio{$pacbio}{'Next'}{'RefStrand'};
                     $ctg_ctg_chain{$right_ctg}{$left_ctg}{2}=$CtgPairPacbio{$pacbio}{'First'}{'RefStrand'};

                     $ctg_ctg_line{$left_ctg}{$right_ctg}=$left_ctg."-".$right_ctg."-".$pacbio;
                     $ctg_ctg_line{$right_ctg}{$left_ctg}=$right_ctg."-".$left_ctg."-".$pacbio;
                     $best_pacbio=$pacbio;
             }
        }
        my $First_start=0;
        my $First_end=0;
        my $Next_start=0;
        my $Next_end=0;
        print ORI "$left_ctg\t$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}\t$ctg_ctg_chain{$left_ctg}{$right_ctg}{1}\t$ctg_ctg_ori{$left_ctg}{$right_ctg}{1}\t$right_ctg\t$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}\t$ctg_ctg_chain{$left_ctg}{$right_ctg}{2}\t$ctg_ctg_ori{$left_ctg}{$right_ctg}{2}\t$ctg_ctg_count{$left_ctg}{$right_ctg}\n";
        
        if($CtgPairPacbio{$best_pacbio}{'First'}{'RefStrand'}==0){
             $First_start=$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'};
             $First_end=$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'};
        }
        elsif($CtgPairPacbio{$best_pacbio}{'First'}{'RefStrand'}==1){
             $First_start=$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'};
             $First_end=$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'};
        }
        if($CtgPairPacbio{$best_pacbio}{'Next'}{'RefStrand'}==0){
             $Next_start=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'};
             $Next_end=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'};
        }
        elsif($CtgPairPacbio{$best_pacbio}{'Next'}{'RefStrand'}==1){
             $Next_start=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'};
             $Next_end=$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'};
        }
        print PAIR "$left_ctg\t$CtgPairPacbio{$best_pacbio}{'First'}{'RefLength'}\t$right_ctg\t$CtgPairPacbio{$best_pacbio}{'Next'}{'RefLength'}\t$ctg_ctg_count{$left_ctg}{$right_ctg}\t$best_pacbio-$ctg_ctg_chain{$left_ctg}{$right_ctg}{1}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefStart'}-$CtgPairPacbio{$best_pacbio}{'First'}{'RefEnd'}-$ctg_ctg_chain{$left_ctg}{$right_ctg}{2}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefStart'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'RefEnd'}\t$best_pacbio-$CtgPairPacbio{$best_pacbio}{'First'}{'QryLength'}-$CtgPairPacbio{$best_pacbio}{'First'}{'QryStart'}-$CtgPairPacbio{$best_pacbio}{'First'}{'QryEnd'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'QryStart'}-$CtgPairPacbio{$best_pacbio}{'Next'}{'QryEnd'}\n";        
    }
    return(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line);
    close(ORI);
    close(PAIR);
}

sub initializeGraph{

    my $ContigNameListFile = shift;
    my $ContigGraph = Graph::Undirected->new;
    my $EasyGraph   = Graph::Easy->new( undirected => 1 );

    open(IN, $ContigNameListFile ) or die("Cannot open $ContigNameListFile\n");

    while ( my $line = <IN> ) {

        chomp $line;

        # contig0001_size5554
        my $Contig = $line;
#        $ContigGraph -> add_vertex( $Contig );
#        $EasyGraph   -> add_node( $Contig );
        
    }      
    return $ContigGraph, $EasyGraph ;
}
sub drawEasyGraph{

    my $EasyGraph        = shift;
    my $MergedContigFile = shift;
    my $SubNum           = shift;
    my $graphviz = $EasyGraph -> as_graphviz();
#    my $BACName = ( $MergedContigFile =~ m/TUG.(\S+).pacmerge(\d+)./ )[0];
#    my $PacmergeNum = ( $MergedContigFile =~ m/TUG.(\S+).pacmerge(\d+)./ )[1];
    open my $DOT, "|dot -Tpng -o $MergedContigFile.png" or die ("Cannot open pipe to dot: $!");
    print $DOT $graphviz;
    close $DOT;

}
sub GenerateRecordHash{

    my $inputfile = shift;
    my %NodeHash = ();
    my %EdgeHash = ();
    my %NodePairsHash = ();
    open(IN, $inputfile ) or die("Cannot open $inputfile\n");
    while ( my $line = <IN> ) {
              chomp $line;
              my @Columns   = split( /\s+/, $line );
              my $RefName   = $Columns[4];
              my $QryName   = $Columns[5];
              my $Score     = $Columns[3];

              if ( ! exists ( $NodeHash{$RefName}{$QryName} ) ) {
                     for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                             my $ColumnName = $ColumnNameArray[$Idx];
                             $NodeHash{$RefName}{$QryName}{$ColumnName} = $Columns[$Idx];
                     }
              }else{
                      if ( $Score < $NodeHash{$RefName}{$QryName}{'Score'} ) {

                              for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                                      my $ColumnName = $ColumnNameArray[$Idx];
                                      $NodeHash{$RefName}{$QryName}{$ColumnName} = $Columns[$Idx];
                              }
                       }
              } 
              if ( ! exists ( $EdgeHash{$QryName}{$RefName} ) ) {
                       for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                            my $ColumnName = $ColumnNameArray[$Idx];
                            $EdgeHash{$QryName}{$RefName}{$ColumnName} = $Columns[$Idx];
                       }
              }else{
                        if ( $Score < $EdgeHash{$QryName}{$RefName}{'Score'} ) {
                            for ( my $Idx = 0 ; $Idx < scalar( @Columns ) ; $Idx ++ ) {
                                 my $ColumnName = $ColumnNameArray[$Idx];
                                 $EdgeHash{$QryName}{$RefName}{$ColumnName} = $Columns[$Idx];
                            }
                        }
              }

      }
      close IN;

      return ( \%NodeHash, \%EdgeHash );
}

sub GenerateNodePairHash{

    my $EdgeHashRef = shift;
    my %EdgeHash = %{ $EdgeHashRef };
    my %NodePairsHash = ();
    foreach my $QryName ( keys %EdgeHash ) {
         if ( ( scalar( keys %{ $EdgeHash{$QryName} } ) > 1 ) ) {
            my %QryRefNameHash     = %{ $EdgeHash{$QryName} };
            my @QryStartHashArray  = ( sort { $QryRefNameHash{$a}{'QryStart'} <=> $QryRefNameHash{$b}{'QryStart'} } keys %QryRefNameHash );
           
            for ( my $QryIdxFst = 0; $QryIdxFst < scalar( @QryStartHashArray ) ; $QryIdxFst ++ ) {
                for ( my $QryIdxNxt = $QryIdxFst + 1 ; $QryIdxNxt < scalar( @QryStartHashArray ) ; $QryIdxNxt ++ ) {

                    my $QryStartScoreFirst = $QryStartHashArray[$QryIdxFst];
                    my $QryStartScoreNext  = $QryStartHashArray[$QryIdxNxt];

                    my %RecordFirstHash = %{ $QryRefNameHash{$QryStartScoreFirst} } ;
                    my %RecordNextHash  = %{ $QryRefNameHash{$QryStartScoreNext}  };

                    if ( $RecordFirstHash{'QryPos'} eq 'head' && $RecordNextHash{'QryPos'} eq 'head' ) {
                    }elsif ( $RecordFirstHash{'QryPos'} eq 'tail' && $RecordNextHash{'QryPos'} eq 'tail' ) {
                    }else{
                        my ( $ContigNameFirst, $ContigNameNext ) = sort ( $RecordFirstHash{'RefName'} , $RecordNextHash {'RefName'} );
                        
                        
                      
                        for ( my $Idx = 0 ; $Idx < scalar( @ColumnNameArray ) ; $Idx ++ ) {
                            my $ColumnName = $ColumnNameArray[$Idx];
                            if($ContigNameFirst eq $RecordFirstHash{'RefName'}){
                                $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'} {$ColumnName}  = $RecordFirstHash{$ColumnName};
                                $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}  {$ColumnName}  = $RecordNextHash{$ColumnName};
                            }
                            elsif($ContigNameFirst ne $RecordFirstHash{'RefName'}){
                                $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'} {$ColumnName}  = $RecordNextHash{$ColumnName};
                                $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}  {$ColumnName}  = $RecordFirstHash{$ColumnName};
                            }

                            if($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefPos'} eq 'head' ){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefPos'}="left";  
                            }
                            elsif($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefPos'}eq 'tail'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefPos'}="right";
                            }
                            if($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefPos'} eq 'head'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefPos'}="left";
                            }
                            elsif($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefPos'} eq 'tail'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefPos'}="right";
                            }
                            if($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefStrand'} eq 'reverse'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefStrand'}=1;
                            }
                            elsif($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefStrand'} eq 'forward'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'First'}{'RefStrand'}=0;
                            }
                            if($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefStrand'} eq 'reverse'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefStrand'}=1;
                            }
                            elsif($NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefStrand'} eq 'forward'){
                                   $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {$QryName} {'Next'}{'RefStrand'}=0;
                            }
                        }
                        

#                        $NodePairsHash{ $ContigNameFirst."-".$ContigNameNext } {'EdgeNum'} ++ ;

#                        $ContigGraph -> add_edge ( $ContigNameFirst, $ContigNameNext );
#                        $EasyGraph   -> add_edge ( $ContigNameFirst, $ContigNameNext );
                     }
                  }
          }
          }else {
          }
    }
    return ( \%NodePairsHash );
}


#&Finding_best_pathway(\%ctg_ctg_count, \%ctg_ctg_overlap_len, \%ctg_ctg_ori, \%ctg_ctg_chain, \%ctg_ctg_line, \%NodePairsHash );
sub Finding_best_pathway{
    my $ctg_pair_count = shift;
    my $ctg_ctg_overlap = shift;
    my $ctg_ori = shift;
    my $ctg_chain = shift;
    my $ctg_line =shift ;
    my $PairsHash = shift;
    my %ctg_ctg_count = %{$ctg_pair_count};
    my %ctg_ctg_overlap = %{$ctg_ctg_overlap};
    my %ctg_ctg_ori = %{$ctg_ori};
    my %ctg_ctg_chain = %{$ctg_chain};
    my %ctg_ctg_line = %{$ctg_line};
    my %NodePairsHash = %{$PairsHash};
    my @cluster_line=();                 #ctg clusters  
    open OUT,">ctg_clusters.txt" or die $!;
    #首先判断一对contig中，是否有一个contig是最好的，如果两个contig都不是最好的，过；否则延伸

    my %existed_pairs=();
    my %all_existed_ctgs=();
    my %all_existed_pairs=();
    foreach my $left_ctg (keys %ctg_ctg_line){
          my $temp_ctg=$ctg_ctg_line{$left_ctg};
          foreach my $right_ctg (keys %$temp_ctg){
             next if(exists $existed_pairs{$left_ctg}{$right_ctg});
             my $left_ctg_order=1;
             my $right_ctg_order=1;
             ################ 判断左边contig是否为其同方向中最好的 ##############################################
             foreach my $left_key1 (keys %ctg_ctg_line){
                 my $left_temp_ctg=$ctg_ctg_line{$left_key1};
                 foreach my $left_key2 (keys %$left_temp_ctg){
                     next if($left_ctg eq $left_key1 && $right_ctg eq $left_key2);
                     if($left_ctg eq $left_key2){
                         if($ctg_ctg_ori{$left_key1}{$left_key2}{2} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}){
                              if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$left_key1}{$left_key2}){
                                  $left_ctg_order++;
                              }
                              elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$left_key1}{$left_key2}){
                                 if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$left_key1}{$left_key2}){
                                     $left_ctg_order++;
                                 }
                              }
                         }
                     }
                 }
              }
              my $left_add=$ctg_ctg_line{$left_ctg};
              foreach my $left_key (keys %$left_add){
                    next if($left_key eq $right_ctg);
                    if($ctg_ctg_ori{$left_ctg}{$left_key}{1} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{1}){
                       if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$left_ctg}{$left_key}){
                           $left_ctg_order++;
                       }
                       elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$left_ctg}{$left_key}){
                           if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$left_ctg}{$left_key}){
                               $left_ctg_order++;
                           }
                       }
                    }
               }
               ##############判断右边contig是否为其同方向中最好的contig #########################################
               foreach my $right_key1 (keys %ctg_ctg_line){
                     my $right_temp_ctg=$ctg_ctg_line{$right_key1};
                     foreach my $right_key2 (keys %$right_temp_ctg){
                         next if($left_ctg eq $right_key1 && $right_ctg eq $right_key2);
                         if($right_ctg eq $right_key2){
                             if($ctg_ctg_ori{$right_key1}{$right_key2}{2} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}){
                                 if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$right_key1}{$right_key2}){
                                    $right_ctg_order++;
                                 }
                                 elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$right_key1}{$right_key2}){
                                     if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$right_key1}{$right_key2}){
                                         $right_ctg_order++;
                                     }
                                 }
                              }
                          }
                      }
                }
                my $right_add=$ctg_ctg_line{$right_ctg};
                foreach my $right_key (keys %$right_add){
                      next if($right_key eq $left_ctg);
                      if($ctg_ctg_ori{$right_ctg}{$right_key}{1} eq $ctg_ctg_ori{$left_ctg}{$right_ctg}{2}){
                           if($ctg_ctg_count{$left_ctg}{$right_ctg}<$ctg_ctg_count{$right_ctg}{$right_key}){
                                  $right_ctg_order++;
                           }
                           elsif($ctg_ctg_count{$left_ctg}{$right_ctg}==$ctg_ctg_count{$right_ctg}{$right_key}){
                              if($ctg_ctg_overlap{$left_ctg}{$right_ctg}<$ctg_ctg_overlap{$right_ctg}{$right_key}){
                                  $right_ctg_order++;
                              }
                           }
                       }
                }
                next if($right_ctg_order>1 || $left_ctg_order>1);
                next if(exists $all_existed_pairs{$left_ctg}{$right_ctg});
                next if(exists $all_existed_ctgs{$left_ctg} || exists $all_existed_ctgs{$right_ctg});
                $all_existed_pairs{$left_ctg}{$right_ctg}=0;
                $all_existed_ctgs{$left_ctg}=0;
                $all_existed_ctgs{$right_ctg}=0;
#                print "$left_ctg\t$right_ctg\n";
                #################### 结束判断 ################################################################
         
                ################### 以上面找到的一对contig为起始 分别左右延伸寻找最好的contig ################ 
                my @extend_line=();
                push(@extend_line,$left_ctg);
                push(@extend_line,$right_ctg);
                my $left_extend=$left_ctg;
                my $right_extend=$right_ctg;

                my %existed_ctg=();
                $existed_ctg{$left_extend}=0;
                $existed_ctg{$right_extend}=0;
         
                my %extended_ctg=();
                my %existed_ctg_ori=();
                $existed_ctg_ori{$left_ctg}=$ctg_ctg_ori{$left_ctg}{$right_ctg}{1};
                $existed_ctg_ori{$right_ctg}=$ctg_ctg_ori{$left_ctg}{$right_ctg}{2};         
                my $left_extend_ctg="1";
                my $left_extend_count=0;
                my $left_overlap=0;
                while($left_extend_ctg ne ""){
                    $left_extend_ctg="";
                    $left_extend_count=0;
                    $left_overlap=0;
                    $left_extend=$extend_line[0];
#                    last if(exists $extended_ctg{$left_extend});      #判断该contig是否已经延伸过
                    $extended_ctg{$left_extend}=0;
                    foreach my $key1 (keys %ctg_ctg_line){
                          my $temp=$ctg_ctg_line{$key1};
                          foreach my $key2 (keys %$temp){
                              if($key2 eq $left_extend){
                                  if($ctg_ctg_ori{$key1}{$key2}{2} ne $ctg_ctg_ori{$left_extend}{$extend_line[1]}{1}){
                                      if($ctg_ctg_count{$key1}{$key2}>$left_extend_count){
                                          $left_extend_ctg=$key1;
                                          $left_extend_count=$ctg_ctg_count{$key1}{$key2};
                                          $left_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                      }
                                      elsif($ctg_ctg_count{$key1}{$key2}==$left_extend_count){
                                          if($ctg_ctg_overlap{$key1}{$key2}>$left_overlap){
                                               $left_extend_ctg=$key1;
                                               $left_extend_count=$ctg_ctg_count{$key1}{$key2};
                                               $left_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                           }
                                      }
                                 }
                              }
                           }
                     }
                     my $temp=$ctg_ctg_line{$left_extend};
                     foreach my $key10 (keys %$temp ){
                         if($ctg_ctg_ori{$left_extend}{$key10}{1} ne $ctg_ctg_ori{$left_extend}{$extend_line[1]}{1}){
                             if($ctg_ctg_count{$left_extend}{$key10}>$left_extend_count){
                                  $left_extend_ctg=$key10;
                                  $left_extend_count=$ctg_ctg_count{$left_extend}{$key10};
                                  $left_overlap=$ctg_ctg_overlap{$left_extend}{$key10};
                             }
                             elsif($ctg_ctg_count{$left_extend}{$key10}==$left_extend_count){
                                  if($ctg_ctg_overlap{$left_extend}{$key10}>$left_overlap){
                                     $left_extend_ctg=$key10;
                                     $left_extend_count=$ctg_ctg_count{$left_extend}{$key10};
                                     $left_overlap=$ctg_ctg_overlap{$left_extend}{$key10};
                                  }
                              }
                          }
                      }
                      if(exists $all_existed_ctgs{$left_extend_ctg}){
                          last;
                      }
                      if(exists $existed_ctg{$left_extend_ctg}){
                          $left_extend_ctg="";
                      }
                      elsif(!exists $existed_ctg{$left_extend_ctg} && $left_extend_ctg ne ""){
                          unshift(@extend_line,$left_extend_ctg);
                          $all_existed_ctgs{$left_extend_ctg}=0;
                          $all_existed_pairs{$left_extend_ctg}{$left_extend}=0;
                          $existed_pairs{$left_extend_ctg}{$left_extend}=0;
                          $existed_pairs{$left_extend_ctg}{$left_extend}=0;
                      }
                      $existed_ctg{$left_extend_ctg}=0;                
                }
        

                ###################### 开始延伸右边 #############################################################
                my $right_extend_ctg="1";
                my $right_extend_count=0;
                my $right_overlap=0;  
                while($right_extend_ctg ne ""){
                    $right_extend_ctg="";
                    $right_extend_count=0;
                    $right_overlap=0;
                    $right_extend=$extend_line[-1];
#                    last if(exists $extended_ctg{$right_extend});          
                    $extended_ctg{$right_extend}=0;
                    foreach my $key1 (keys %ctg_ctg_line){
                        my $temp=$ctg_ctg_line{$key1};
                        foreach my $key2 (keys %$temp){
                           if($key2 eq $right_extend){
                               if($ctg_ctg_ori{$key1}{$key2}{2} ne $ctg_ctg_ori{$right_extend}{$extend_line[-2]}{1}){
                                    if($ctg_ctg_count{$key1}{$key2}>$right_extend_count){
                                               $right_extend_ctg=$key1;
                                               $right_extend_count=$ctg_ctg_count{$key1}{$key2};
                                               $right_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                     }
                                     elsif($ctg_ctg_count{$key1}{$key2}==$right_extend_count){
                                               if($ctg_ctg_overlap{$key1}{$key2}>$right_overlap){
                                                  $right_extend_ctg=$key1;
                                                  $right_extend_count=$ctg_ctg_count{$key1}{$key2};
                                                  $right_overlap=$ctg_ctg_overlap{$key1}{$key2};
                                               }
                                      }
                                }
                            }
                         }
                    }
                    my $temp=$ctg_ctg_line{$right_extend};
                    foreach my $key1 (keys %$temp ){
                        next if($key1 eq $extend_line[-2]);
                        if($ctg_ctg_ori{$right_extend}{$key1}{1} ne $ctg_ctg_ori{$right_extend}{$extend_line[-2]}{1}){
                            if($ctg_ctg_count{$right_extend}{$key1}>$right_extend_count){
                                 $right_extend_ctg=$key1;
                                 $right_extend_count=$ctg_ctg_count{$right_extend}{$key1};
                                 $right_overlap=$ctg_ctg_overlap{$right_extend}{$key1};
                            }
                            elsif($ctg_ctg_count{$right_extend}{$key1}==$right_extend_count){
                                 if($ctg_ctg_overlap{$right_extend}{$key1}>$right_overlap){
                                     $right_extend_ctg=$key1;
                                     $right_extend_count=$ctg_ctg_count{$right_extend}{$key1};
                                     $right_overlap=$ctg_ctg_overlap{$right_extend}{$key1};
                                 }
                            }
                         }
                     }
                     if(exists $all_existed_ctgs{$right_extend_ctg}){
                          last;
                      }
                     if(exists $existed_ctg{$right_extend_ctg}){
                         $right_extend_ctg="";
                     }
                     elsif(!exists $existed_ctg{$right_extend_ctg} && $right_extend_ctg ne ""){
                         push(@extend_line,$right_extend_ctg);
                         $all_existed_ctgs{$right_extend_ctg}=0;
                         $all_existed_pairs{$right_extend}{$right_extend_ctg}=0;
                         $existed_pairs{$right_extend}{$right_extend_ctg}=0;
                         $existed_pairs{$right_extend_ctg}{$right_extend}=0;
                     }
                     $existed_ctg{$right_extend_ctg}=0;
                }
                my $ctg_cluster=$extend_line[0];

                for(my $i=1;$i<@extend_line;$i++){
                      $ctg_cluster=$ctg_cluster."-".$extend_line[$i];
                }
#                print "$left_ctg $right_ctg----$ctg_cluster\n";
                my $sign=0;
                for(my $i=0;$i<@cluster_line;$i++){
                      if($cluster_line[$i]=~/$ctg_cluster/){
                             $sign=1;
                      }
                      elsif($ctg_cluster=~/$cluster_line[$i]/){
                             $cluster_line[$i]=$ctg_cluster;
                             $sign=1;
                      }
                      else{
                             my @content=split "-",$ctg_cluster;
                             my $reverse=join ("-",reverse(@content));
                             if($reverse=~/$cluster_line[$i]/){
                                 $cluster_line[$i]=$reverse;
                                 $sign=1;
                             }
                             elsif($cluster_line[$i]=~/$reverse/){
                                 $sign=1;
                             }
                      }
                 }
                 ################# 检查所有已经找到的通路中是否有包含当前通路中所有contig的通路 ##############
 
                 for(my $j=0;$j<@cluster_line;$j++){
                      my @target_ctg=split "-",$cluster_line[$j];
                      my %ctg_target=();
                      my %ctg_query=();
                      for(my $target=0;$target<@target_ctg;$target++){
                           $ctg_target{$target_ctg[$target]}=0;
                      }
                      for(my $query=0;$query<@extend_line;$query++){
                           $ctg_query{$extend_line[$query]}=0;
                      }
                      my $target_sign=0;
                      foreach my $target_key(keys %ctg_target){
                           if(!exists $ctg_query{$target_key}){
                                    $target_sign=1;
                           }
                      }
                      if($target_sign==0){
                           $cluster_line[$j]=$ctg_cluster;
                           $sign=1;
                      }
                      my $query_sign=0;
                      foreach my $query_key(keys %ctg_query){
                           if(!exists $ctg_target{$query_key}){
                                    $query_sign=1;
                           }
                      }
                      if($query_sign==0){
                           $sign=1;
                      }
                 }
                 ############################ end #######################################################  
                 if($sign==0){
                      push(@cluster_line,$ctg_cluster);
                 }           
          }
    }
#    print "@cluster_line\n";
    for(my $i=0;$i<@cluster_line;$i++){
        my @content=split "-",$cluster_line[$i];
        for(my $j=0;$j<@content;$j++){
            print OUT "$content[$j] ";
        } 
        print OUT "\n";
    }
    close(OUT);
    return (@cluster_line);
}                                          
sub draw_ctg_clusters{
    my $info=shift;
    my @ctg_cluster=@{$info};
    my $ClusterGraph   = Graph::Easy->new( undirected => 1 );
    my %existed_ctg=();
#    print "@ctg_cluster\n";
    open CLUSTER,">cluster_ori.txt" or die $!;
    for(my $i=0;$i<@ctg_cluster;$i++){
        my @content=split "-",$ctg_cluster[$i];
=pod
        for(my $j=0;$j<@content;$j++){
             next if(exists $existed_ctg{$content[$j]});
             $EasyGraph   -> add_node( $content[$j] );
        }
=cut
        for(my $j=0;$j<@content-1;$j++){
             $ClusterGraph -> add_edge ($content[$j],$content[$j+1]);
             my @infor=split "-",$ctg_ctg_line{$content[$j]}{$content[$j+1]};
             print CLUSTER "$content[$j] $ctg_ctg_chain{$content[$j]}{$content[$j+1]}{1} $ctg_ctg_ori{$content[$j]}{$content[$j+1]}{1} $ctg_ctg_chain{$content[$j]}{$content[$j+1]}{2} $ctg_ctg_ori{$content[$j]}{$content[$j+1]}{2} ";
        }
        print CLUSTER "$content[-1]\n";
     }
     &drawEasyGraph ( $ClusterGraph , "ctg_cluster" , '01' );
     
}     
                         
 
open OUT1,">pairs.txt" or die $!;
foreach my $key1 (keys %NodePairsHash){
       my $temp=$NodePairsHash{$key1};
       my @content=split "-",$key1;
       foreach my $pacbio (keys %$temp){
            print OUT1 "$content[0] $NodePairsHash{$key1}{$pacbio}{'First'}{'RefPos'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefStrand'} $NodePairsHash{$key1}{$pacbio}{'First'}{'Identity'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefStart'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefEnd'} $NodePairsHash{$key1}{$pacbio}{'First'}{'RefLength'}-----$content[1] $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefPos'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefStrand'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'Identity'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefStart'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefEnd'} $NodePairsHash{$key1}{$pacbio}{'Next'}{'RefLength'}\n";
       }
}

