#!/usr/bin/env perl

# 寻找两端100%比对上1000bp的序列和完全包含的序列，简化节点数
# 做连接和删除后重新比对
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Perl;

my $BlasrAlignFile       = shift;
my $ContigFile           = shift;
my $MergedContigFile     = shift;
my $MinOverhangLength    = 500;
my $MinOverlapLength     = shift;

my %DisposeHash = ();
my %MergeHash   = ();

open(IN, $BlasrAlignFile ) or die("Cannot open $BlasrAlignFile\n");

while ( my $line = <IN> ) {

	chomp $line; 
    
    next if ( $line =~ m/nCells/ );
    
    # qName	               tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
    # contig0015_size25589 contig0003_size29117 0       0       -110660 100                 0       22132   29117   3457    25589   25589 464738
    # contig0030_size1680  contig0032_size1680  0       0       -8400   100                 0       1680    1680    0       1680    1680  35246
    my @Columns   = split( /\s+/, $line );
   
    my $score     = $Columns[4];
    my $percentSimilarity   = $Columns[5];
    
    # Illumina Contig
    my $tName     = $Columns[1];
    my $tStrand   = $Columns[3];
    my $tStart    = $Columns[6];
    my $tEnd      = $Columns[7];
    my $tLength   = $Columns[8];
    # Illumina Contig
    my $qName     = $Columns[0];
    my $qStrand   = $Columns[2];
    my $qStart    = $Columns[9];
    my $qEnd      = $Columns[10];
    my $qLength   = $Columns[11];
    
    # 比对结果1K以下的之后处理
    
    my $RefLeftOverhang  = 0;
    my $RefRightOverhang = 0;

    my $QryLeftOverhang  = 0 ;
    my $QryRightOverhang = 0 ;

    $RefLeftOverhang  = $tStart;
    $RefRightOverhang = ( $tLength - $tEnd );
    $QryLeftOverhang  = $qStart;
    $QryRightOverhang = ( $qLength - $qEnd );

    # 包含关系检测
    if ( $RefLeftOverhang <=$MinOverhangLength && $RefRightOverhang <= $MinOverhangLength &&
         $QryLeftOverhang <=$MinOverhangLength && $QryRightOverhang <= $MinOverhangLength    ) {
        # 完全相等，扔一边
        print $line."\n\n";
        $DisposeHash{$tName} = 0;
        next;
    }elsif ( $RefLeftOverhang <=$MinOverhangLength && $RefRightOverhang <=$MinOverhangLength ) {
        # Ref包含,舍弃Ref
        print $line."\n";
        $DisposeHash{$tName} = 0;
        next;
    }elsif( $QryLeftOverhang <=$MinOverhangLength && $QryRightOverhang <=$MinOverhangLength ) {
        # Qry包含,舍弃Qry
        print $line."\n";
        $DisposeHash{$qName} = 0;
        next;
    }

    # 比对关系检测，寻找两端100%比对上1000bp的序列
    if ( $percentSimilarity >=98 && abs( $qEnd - $qStart ) >= $MinOverlapLength ) {
    
        if ( $qStrand == 0 && $ tStrand == 0 ) {

            # 相同方向，head
            #            ======================>
            #            |||||||||||||||
            #       ------------------->
            # qName	                tName	                qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
            # contig0015_size25589  contig0003_size29117    0       0       -110660 100                 0       22132   29117   3457    25589   25589   464738
            if ( ( $qEnd == $qLength ) && ( $qStart  > 0 ) &&
                 ( $tStart == 0 )      && ( $tEnd    < $tLength  ) ) {
                print $line."\n";
                print "$qName\tf\t0\t$qEnd\t$tName\tf\t$tEnd\t$tLength\n";
                $MergeHash{abs($score)} = "$qName\tf\t0\t$qEnd\t$tName\tf\t$tEnd\t$tLength";
            # 相同方向，tail
            #         ==========================>
            #                     |||||||||||||||
            #                     --------------------->
            # qName	            tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
            # 
            }elsif ( ( $qEnd   < $qLength ) && ( $qStart  == 0 ) &&
                     ( $tStart > 0 )        && ( $tEnd    == $tLength ) ) {
                print $line."\n";
                print "$tName\tf\t0\t$tEnd\t$qName\tf\t$qEnd\t$qLength\n";
                $MergeHash{abs($score)} = "$tName\tf\t0\t$tEnd\t$qName\tf\t$qEnd\t$qLength";
            }

        }elsif ( $qStrand == 0 && $ tStrand == 1 ) {

            $RefLeftOverhang  = $tStart;
            $RefRightOverhang = ( $tLength - $tEnd );
            $QryLeftOverhang  = $qStart;
            $QryRightOverhang = ( $qLength - $qEnd );

            # 不同方向，head
            #    <======================
            #            |||||||||||||||          <-----
            #            -------------------->
            #　qName	            tName	                qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
            #　contig0006_size2382  contig0018_size2968      0       1       -6681  98.9759             1601    2968    2968    0       1367    2382    28711
            if ( ( $qEnd   < $qLength ) && ( $qStart  == 0 ) &&
                 ( $tStart > 0 )        && ( $tEnd    == $tLength  ) ) {
                print $line."\n";
                print "$tName\tr\t0\t$tLength\t$qName\tf\t$qEnd\t$qLength\n";
                $MergeHash{abs($score)} = "$tName\tr\t0\t$tLength\t$qName\tf\t$qEnd\t$qLength";
                
            # 不同方向，tail
            #          <======================
            #          |||||||||||||||            <-----
            #    -------------------->
            # qName	                tName	            qStrand	tStrand	score	percentSimilarity	tStart	tEnd	tLength	qStart	qEnd	qLength	nCells
            # contig0023_size1033   contig0009_size480  0       0       -1490   100                 0       298     480     735     1033    1033    6224
            }elsif ( ( $qEnd   == $qLength ) && ( $qStart  > 0 ) &&
                     ( $tStart == 0 )        && ( $tEnd    < $tLength  ) ) {
                print $line."\n";
                print "$qName\tf\t0\t$qEnd\t$tName\tr\t$tEnd\t$tLength\n";
                $MergeHash{abs($score)} = "$qName\tf\t0\t$qEnd\t$tName\tr\t$tEnd\t$tLength";
            }
            
        }
    }
    
}

close IN;

my %SeqHash = ();

# 读入序列
my $in  = Bio::SeqIO->new(-file => "<$ContigFile" , '-format' => 'Fasta');
            
    while ( my $seq = $in->next_seq() ) {
        
        my $Contig       = $seq->id;
        my $ContigLength = $seq->length;
        my $ContigSeq    = $seq->seq;
                
        if ( ! exists( $DisposeHash{$Contig} ) ) {
            $SeqHash{$Contig} = $ContigSeq;
        }else{
            # skip
        }
        
    }

open(OUT, ">$MergedContigFile" ) or die("Cannot open $MergedContigFile\n");

my %UsedHash = ();

# 合并序列
foreach my $Score ( sort { $b <=> $a } keys %MergeHash ) {
    
    my $MergeRecord = $MergeHash{$Score};
    my @Columns = split( /\t/, $MergeRecord );
    
    my $LeftName   = $Columns[0];
    my $LeftStrand = $Columns[1];
    my $LeftStart  = $Columns[2];
    my $LeftEnd    = $Columns[3];

    my $RightName   = $Columns[4];
    my $RightStrand = $Columns[5];
    my $RightStart  = $Columns[6];
    my $RightEnd    = $Columns[7];
    
    next if ( exists( $DisposeHash{$LeftName} ) || exists( $DisposeHash{$RightName} ) );
    # next if ( exists( $UsedHash{$LeftName."\t".$RightName} ) || exists( $UsedHash{$RightName."\t".$LeftName} ) );
    next if ( ! exists( $SeqHash{$LeftName} ) || ! exists( $SeqHash{$RightName} ) );
    
    my $MergeName = $LeftName."_".$RightName;
    my $MergeSeq  = '';
    
    if ( $LeftStrand eq 'f' ) {
       $MergeSeq .= substr( $SeqHash{$LeftName}, $LeftStart, ( $LeftEnd - $LeftStart ) );
    }else{
       my $revcom = reverse $SeqHash{$LeftName};
       $revcom =~ tr/ACGTacgt/TGCAtgca/;
       $MergeSeq .= substr( $revcom, $LeftStart, ( $LeftEnd - $LeftStart ) );
    }
    delete $SeqHash{$LeftName};
    
    if ( $RightStrand eq 'f' ) {
       $MergeSeq .= substr( $SeqHash{$RightName}, $RightStart, ( $RightEnd - $RightStart ) );
    }else{
       my $revcom = reverse $SeqHash{$RightName};
       $revcom =~ tr/ACGTacgt/TGCAtgca/;
       $MergeSeq .= substr( $revcom, $RightStart, ( $RightEnd - $RightStart ) );
    }
    delete $SeqHash{$RightName};
    
    print OUT ">$MergeName\n";
    print OUT "$MergeSeq\n";
    
    # $UsedHash{$LeftName."\t".$RightName} = 0;
    
}

foreach my $Contig ( sort{ $a cmp $b } keys %SeqHash ) {
    my $ContigSeq = $SeqHash{$Contig};
    print OUT ">$Contig\n";
    print OUT "$ContigSeq\n";
}

close OUT;

print "Disposed\n";
foreach my $DisposeRecord ( keys %DisposeHash ) {
    print "$DisposeRecord\n";
}

# print "Used\n";
# foreach my $UsedRecord ( keys %UsedHash ) {
    # print "$UsedRecord\n";
# }
