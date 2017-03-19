#!/usr/bin/perl
#Author: huilong du 
#Node: merged ctg into scaffold
use warnings;
use strict;

my $infile1=shift;             #list_all_pac.txt
my $infile2=shift;             #list_all_ctg.txt

open IN1,"<$infile1" or die $!;
open IN2,"<$infile2" or die $!;

my %pac_file=();
my %ctg_file=();
my $sign="";

while(<IN1>){
     chomp;
     my $line=$_;
     my @content=split "/",$line;
     $content[-1]=~/(\S+).fasta/;
     my $sign=$1;
     $pac_file{$sign}=$line;
   
}

while(<IN2>){
      chomp;
      my $line=$_;
      my @content=split "/",$line;
      $content[-1]=~/(\S+).ctg.fasta/;
      my $sign=$1;
      $ctg_file{$sign}=$line;
}
close(IN1);
close(IN2);

my $count=0;

foreach my $key (keys %ctg_file){
      next if(!exists $pac_file{$key});
      system("mkdir $key");      
      open OUT,">./$key/$count.pbs" or die $!;
      if($count%4<=2){
          print OUT "#PBS -N $count-$key-wheat-merge_ctg
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q hldu
#PBS -l mem=5gb
cd \$PBS_O_WORKDIR
";
       }
       elsif($count%4==3){
          print OUT "#PBS -N $count-$key-wheat-merge_ctg
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q hldu_vvl
#PBS -l mem=5gb
cd \$PBS_O_WORKDIR
";
       }
          
       print OUT "cp /public/share/hldu/wheat_bac_pacbio/script_pipeline/*.pl ./$key/\n";
       print OUT "cd ./$key/\n";
       print OUT "time blasr \\
                       $ctg_file{$key}\\
                       $ctg_file{$key}\\
                       -minMatch 14 \\
                       -bestn 10 \\
                       -nproc 4\\
                       -header \\
                       -out $key.merge0.txt\n";
            print OUT "sed 's/\\/0_[0-9]\\+//' $key.merge0.txt | awk '{if(\$1!=\$2 && \$6 > 95 ){print}}' >$key.merge0.Filtered.txt\n";
            print OUT "time perl 01-BLASRDeOverlap.pl ./$key.merge0.Filtered.txt $ctg_file{$key} $key.merge1.fasta 1000\n";
            print OUT "time blasr \\
        ./$key.merge1.fasta \\
        ./$key.merge1.fasta \\
        -minMatch 14 \\
        -bestn 10 \\
        -nproc 4 \\
        -header \\
        -out ./$key.merge1.txt
";
             print OUT "sed 's/\\/0_[0-9]\\+//' ./$key.merge1.txt | awk '{if(\$1!=\$2 && \$6 > 95 ){print}}' > ./$key.merge1.Filtered.txt\n";
             print OUT "time perl 01-BLASRDeOverlap.pl $key.merge1.Filtered.txt $key.merge1.fasta $key.merge2.fasta 500\n";
             print OUT "time blasr \\
        ./$key.merge2.fasta \\
        ./$key.merge2.fasta \\
        -minMatch 14 \\
        -bestn 10 \\
        -nproc 4 \\
        -header \\
        -out $key.merge2.txt
";
              print OUT "sed 's/\\/0_[0-9]\\+//' ./$key.merge2.txt | awk '{if(\$1!=\$2 && \$6 > 95 ){print}}' > ./$key.merge2.Filtered.txt\n";
              print OUT "time perl 01-BLASRDeOverlap.pl ./$key.merge2.Filtered.txt $key.merge2.fasta ./$key.merge3.fasta 300\n";
              print OUT "cat ./$key.merge3.fasta| grep \">\" | sed \'s/>//\' > ./$key.merge3.list\n";
              print OUT "sawriter ./$key.merge3.fasta\n";
              print OUT "time blasr \\
        $pac_file{$key} \\
        $key.merge3.fasta \\
        -sa $key.merge3.fasta.sa \\
        -minMatch 14 \\
        -bestn 10 \\
        -nproc 4 \\
        -header \\
        -out $key.pacmerge_blasr.txt
";
              print OUT "time perl 02-PacbioGAPFilter.pl $key.pacmerge_blasr.txt 500 95 > 02-PacbioGAPFilter.txt\n";
              print OUT "time perl 03-PacbioGAPLinker.pl 02-PacbioGAPFilter.txt   500        > 03-PacbioGAPLinker.txt\n";
              print OUT "time perl remove_include.pl 03-PacbioGAPLinker.txt 03-PacbioGAPLinker_filtered.txt\n";

              print OUT "time perl construct_graph_final.pl 03-PacbioGAPLinker_filtered.txt $key.merge3.list $key\n";
              print OUT "perl make_ctg_line.pl cluster_ori.txt cluster_ori_same_chain.txt\n";
              print OUT "perl make_junction_by_pos.pl ctg_pairs.txt ctg_ctg_ori.txt cluster_ori_same_chain.txt cluster_ori_same_chain_pos.txt\n";
              print OUT "perl extract_ctg_infor_for_seq.pl cluster_ori_same_chain_pos.txt cluster_ori_same_chain_pos_for_seq.txt\n";
              print OUT "perl extract_seq_by_pos_test.pl cluster_ori_same_chain_pos_for_seq.txt $key.merge3.fasta $pac_file{$key} $key.scaffold.fasta\n"; 
              print OUT "mv ctg_cluster.png $key-cluster.png\n";
              print OUT "cat $key.scaffold.fasta non-connected-ctg.fasta>$key-bac.fasta\n";
              close(OUT);  
              my $Jobs=`qstat -u hldu | wc -l`;

              while( $Jobs > 450 ){
                       print "JOB Remain $Jobs\n";
                       sleep 30;
                       print `date`;
                       $Jobs=`qstat -u hldu | wc -l`;
              }

              my $commond=`qsub ./$key/$count.pbs`;
              $count++;

}
