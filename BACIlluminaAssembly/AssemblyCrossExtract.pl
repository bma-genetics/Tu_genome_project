use strict;
use warnings;
use List::Util qw( min max );
use List::MoreUtils qw( minmax firstidx first_index lastidx last_index indexes);

my $DEBUG = 0;

my $HorizontalFile = 'Horizontal.txt';
my $VerticalFile   = 'Vertical.txt';

my @SampleArrayHorizontal = ();
my @SampleArrayVertical   = ();

open( HORFILE, "<$HorizontalFile") || die "cannot open file $HorizontalFile \n" ;
    
    while(<HORFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        push ( @SampleArrayHorizontal, $line );
        
    } # end of while

close( HORFILE );

open( VERFILE, "<$VerticalFile") || die "cannot open file $VerticalFile \n" ;

    my $QueryCount = 0;

    while(<VERFILE>)
    {
        my($line) = $_;
        chomp($line);

        push ( @SampleArrayVertical, $line );
        
    } # end of while

close( VERFILE );

my $PreFixQuery = shift; # 'G2S1';

my $ScoreFile    = "./AssemblyCrossScore/$PreFixQuery.txt";
my $ReadNameFile = "./FilteredFasta/$PreFixQuery-Filtered.fasta";

my $ExtractQuota    = shift;
my $MinimumCoverage = shift;
my $Times           = shift;

my @ExtractFileArray = ();
my @ExtractPoolArray = ();

# 判断输入的样本名称是在垂直方向还是水平方向
if ( grep {/$PreFixQuery/} @SampleArrayVertical ) {
    
    # 垂直方向，打开水平方向文件
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {
        # PSL比对文件
        my ($ExtractFile) = "./AssemblyCrossExtract/$PreFixQuery-$SampleHorizontal.txt";
        # print $ExtractFile."\n";
        open my $FileHandler, ">$ExtractFile";
        push ( @ExtractFileArray, $FileHandler );
        push ( @ExtractPoolArray, $SampleHorizontal );
    }
    
}elsif ( grep {/$PreFixQuery/} @SampleArrayHorizontal ) {
    
    # 水平方向，打开垂直方向文件
    foreach my $SampleVertical ( @SampleArrayVertical ) {
        # PSL比对文件
        my ($ExtractFile) = "./AssemblyCrossExtract/$PreFixQuery-$SampleVertical.txt";
        # print $ExtractFile."\n";
        open my $FileHandler, ">$ExtractFile";
        push ( @ExtractFileArray, $FileHandler );
        push ( @ExtractPoolArray, $SampleVertical );
    }
    
}else{
    print STDERR "No Such SampleName : $PreFixQuery\n";
    exit;
}

my %ReadCountHash = ();

# 读取ReadName,建立相应的向量数组
open( READFILE, "<$ReadNameFile") || die "cannot open file $ReadNameFile \n" ;

    while(<READFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        if ( $line =~ />(\S+)/ ) {
            my $ReadName = $1;
            $ReadCountHash{$ReadName} = 0;
        }

    } # end of while

close( READFILE );

my $ReadCount = 0;

my %ExtractNumReadCountHash = ();
my $ZeroDepthReadCount = 0;
my $LowDepthReadCount  = 0;
my $ReapeatReadCount   = 0;
my $AllReadCount       = 0;

# 读入Score文件的Read名字和分数
open(INFILE, "<$ScoreFile") || die "cannot open file $ScoreFile \n" ;

    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
		
        $ReadCount++;

        my ( $ReadName, $Scores ) = split( /\t/, $line, 2 );
        
        my @ScoreArray = split( /\t/, $Scores );
        my @BACKScoreArray = @ScoreArray;
        
        # 全部交叉混池中都出现的重复序列去掉
        if ( !defined(scalar( indexes { $_ == 0 } @ScoreArray )) ) {
            # print "Repeat : ".join( " ", @ScoreArray )."\n";
            $ReapeatReadCount ++;
            $AllReadCount ++;
            next;
        }

        my @SortedScoreArray = sort { $b <=> $a } @ScoreArray;
        
        # 全部交叉混池中都是0的去掉
        if ( $SortedScoreArray[0] == 0 ) {
            # print "ZeroDepth : ".join( " ", @ScoreArray )."\n";
            $ZeroDepthReadCount ++;
            $AllReadCount ++;
            next;
        }
        
        # 最大深度达不到最低深度要求的去掉
        if ( $SortedScoreArray[0] <= $MinimumCoverage ) {
            # print "LowDepth : ".join( " ", @ScoreArray )."\n";
            $LowDepthReadCount ++;
            $AllReadCount ++;
            next;
        }

        if ( $DEBUG == 1 ) {
            print "$ReadName\n";
            print join( " ", @ScoreArray )."\n";
            print join( " ", @SortedScoreArray )."\n";
        }
        
        # Coverage小于1，不取log
        # 排序前清除只比对上一次的数组元素, 其他数值取log2
        # for( my $ScoreArrayIndex = 0; $ScoreArrayIndex < scalar( @ScoreArray ); $ScoreArrayIndex ++ ) {
        #     $ScoreArray[$ScoreArrayIndex] = 1 if ( $ScoreArray[$ScoreArrayIndex] == 0 );
        #     $ScoreArray[$ScoreArrayIndex] = log( $ScoreArray[$ScoreArrayIndex] ) / log(2);
        # }
        # 排序后清除只比对上一次的数组元素, 其他数值取log2
        # for( my $ScoreArrayIndex = 0; $ScoreArrayIndex < scalar( @SortedScoreArray ); $ScoreArrayIndex ++ ) {
        #     $SortedScoreArray[$ScoreArrayIndex] = 1 if ( $SortedScoreArray[$ScoreArrayIndex] == 0 );
        #     $SortedScoreArray[$ScoreArrayIndex] = log( $SortedScoreArray[$ScoreArrayIndex] ) / log(2);
        # }
        
        if ( $DEBUG == 1 ) {
            print join( " ", @ScoreArray )."\n";
            print join( " ", @SortedScoreArray )."\n";
        }
        
        # 排序后的分数向后搜索
        # 到零的时候停止搜索
        # 取所有大于零的值的平均值，计算有多少个混池大于平均值
        # 明显比较多的混池个数小于允许的个数时输出到相应的BAC中去
        
        my $ScoreSum = 0;
        my $ScoreNum = 0;
        
        for( my $SortedScoreIndex = 0; $SortedScoreIndex < scalar( @SortedScoreArray ); $SortedScoreIndex ++ ) {
            
            # 分数不为零时累积求平均值
            if  ( $SortedScoreArray[$SortedScoreIndex] != 0 ) {
            
                $ScoreSum = $ScoreSum + $SortedScoreArray[$SortedScoreIndex];
                $ScoreNum ++;
                next;
            
            }elsif  ( $SortedScoreArray[$SortedScoreIndex] == 0 ) {
                
                # 求平均值
                my $AverageScore = $ScoreSum / $ScoreNum ;
                print "$AverageScore\n" if ( $DEBUG == 1 );
                
                # my @Indices01 = indexes { $_ >= $AverageScore } @SortedScoreArray;
                # print join( " ", @Indices01 )."\n";
                
                # my @Indices = indexes { $_ >= $AverageScore } @ScoreArray; # Coverage不用平均值过滤
                my @Indices = indexes { $_ >= $MinimumCoverage } @ScoreArray;
                print join( " ", @Indices )."\n" if ( $DEBUG == 1 );
                
                my $ExtractNum = scalar( @Indices );
                print "$ExtractNum\n" if ( $DEBUG == 1 );
                $ExtractNumReadCountHash{$ExtractNum} ++;
                $AllReadCount ++;
                
                # 明显比较多的混池个数小于允许的个数输出
                if ( $ExtractNum <= $ExtractQuota) {
                    
                    for ( my $OutputIndex = 0; $OutputIndex < $ExtractNum; $OutputIndex ++ ) {
                        
                        print join( " ", @Indices )."\n" if ( $DEBUG == 1 );
                        print $Indices[$OutputIndex]."\n" if ( $DEBUG == 1 );
                        
                        print { $ExtractFileArray[$Indices[$OutputIndex]] } $ReadName."\t";
                        print { $ExtractFileArray[$Indices[$OutputIndex]] } $PreFixQuery."\t";
                        # print { $ExtractFileArray[$Indices[$OutputIndex]] } $ReadCountHash{$ReadName}."\t";
                        print { $ExtractFileArray[$Indices[$OutputIndex]] } $ExtractPoolArray[$Indices[$OutputIndex]]."\t";
                        print { $ExtractFileArray[$Indices[$OutputIndex]] } $BACKScoreArray[$Indices[$OutputIndex]]."\n";
                        
                    }
                    # print "\n";
                }
                last;
                
            } # SortedScoreArray
            
        } # SortedScoreIndex
        
        # print "\n";

    } # end of while

close(INFILE);

# foreach my $FileHandler ( @ExtractFileArray ) {
  # close $FileHandler;
# }
print "ReadName\tOne\tTwo\tThree\tFour\tFive\tSix\tSeven\tEight\tNine\tTen\tZeroDepth\tAllExtract\n"; # if ( $DEBUG == 1 );
print ("$PreFixQuery\t");
print ("$ExtractNumReadCountHash{1}\t");
print ("$ExtractNumReadCountHash{2}\t");
print ("$ExtractNumReadCountHash{3}\t");
print ("$ExtractNumReadCountHash{4}\t");
print ("$ExtractNumReadCountHash{5}\t");
print ("$ExtractNumReadCountHash{6}\t");
print ("$ExtractNumReadCountHash{7}\t");
print ("$ExtractNumReadCountHash{8}\t");
print ("$ExtractNumReadCountHash{9}\t");
print ("$ExtractNumReadCountHash{10}\t");
print ("$ZeroDepthReadCount\t");
print ("$AllReadCount\n");

print ("\n\n");

print "ReadName\tOne\tTwo\tThree\tFour\tFive\tSix\tSeven\tEight\tNine\tTen\tZeroDepth\tAllExtract\n"; # if ( $DEBUG == 1 );
print ("$PreFixQuery\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{1} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{2} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{3} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{4} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{5} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{6} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{7} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{8} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{9} / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ExtractNumReadCountHash{10}/ $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ZeroDepthReadCount         / $AllReadCount ) ); print ("\t");
print ("$AllReadCount\n");

print ("\n\n");

print "LowDepth\tRepeat\tAllReads\n"; # if ( $DEBUG == 1 );
print ("$LowDepthReadCount\t");
print ("$ReapeatReadCount\t");
print ("$ReadCount\n");

print ("\n\n");

print "LowDepth\tRepeat\tAllReads\n"; # if ( $DEBUG == 1 );
printf ( "%.2f", ( $LowDepthReadCount  / $AllReadCount ) ); print ("\t");
printf ( "%.2f", ( $ReapeatReadCount   / $AllReadCount ) ); print ("\t");
print ("$ReadCount\n");

# foreach my $ExtractNum ( sort { $a <=> $b } keys %ExtractNumReadCountHash ) {
    # print STDERR ("$ExtractNum:$ExtractNumReadCountHash{$ExtractNum}\t");
# }

# print join( " ", @ScoreArray )."\n";
# print join( " ", @SortedScoreArray )."\n";
                
# printf "FirstIndex0 %i\n", firstidx { $_ == $SortedScoreArray[0] } @ScoreArray;
# printf "FirstIndex1 %i\n", firstidx { $_ == $SortedScoreArray[1] } @ScoreArray;
# printf "FirstIndex2 %i\n", firstidx { $_ == $SortedScoreArray[2] } @ScoreArray;

# printf "LastIndex0 %i\n", lastidx { $_ == $SortedScoreArray[0] } @ScoreArray;
# printf "LastIndex1 %i\n", lastidx { $_ == $SortedScoreArray[1] } @ScoreArray;
# printf "LastIndex2 %i\n", lastidx { $_ == $SortedScoreArray[2] } @ScoreArray;

# my @x = indexes { $_ == $SortedScoreArray[0] } @ScoreArray;
# print join( " ", @x )."\n";
# my @y = indexes { $_ == $SortedScoreArray[1] } @ScoreArray;
# print join( " ", @y )."\n";
# my @z = indexes { $_ == $SortedScoreArray[2] } @ScoreArray;
# print join( " ", @z )."\n";
