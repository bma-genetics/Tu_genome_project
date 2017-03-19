use strict;
use warnings;
use Bio::SeqIO;

# 计算一个混池中的Reads在交叉混池中出现的覆盖度

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

my $ReadNameFile = "./FilteredFasta/$PreFixQuery-Filtered.fasta";
my $FastAOutFile = "./AssemblyCrossScore/$PreFixQuery.txt";

my @PSLFileArray = ();

# 判断输入的样本名称是在垂直方向还是水平方向
if ( grep {/$PreFixQuery/} @SampleArrayVertical ) {
    
    # 垂直方向，打开水平方向文件
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {
        # PSL比对文件
        my ($PSLFile) = "./AssemblyCrossCount/$PreFixQuery-$SampleHorizontal.txt";
        # print $PSLFile."\n";
        open my $FileHandler, "<$PSLFile";
        push ( @PSLFileArray, $FileHandler );

    }
    
}elsif ( grep {/$PreFixQuery/} @SampleArrayHorizontal ) {
    
    # 水平方向，打开垂直方向文件
    foreach my $SampleVertical ( @SampleArrayVertical ) {
        # PSL比对文件
        my ($PSLFile) = "./AssemblyCrossCount/$PreFixQuery-$SampleVertical.txt";
        # print $PSLFile."\n";
        open my $FileHandler, "<$PSLFile";
        push ( @PSLFileArray, $FileHandler );
   
    }
    
}else{
    print STDERR "No Such SampleName : $PreFixQuery\n";
    exit;
}

my %ReadScoreHash = ();

# 读取ReadName,建立相应的向量数组
open( READFILE, "<$ReadNameFile") || die "cannot open file $ReadNameFile \n" ;

    while(<READFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        if ( $line =~ />(\S+)/ ) {
            my $ReadName = $1;
            my @TempArray = (0) x scalar( @PSLFileArray );
            $ReadScoreHash{$ReadName} = \@TempArray;
        }

        
    } # end of while

close( READFILE );

my $PSLFileCount = 0;

foreach my $FileHandler ( @PSLFileArray ) {
    
    # 读取Count结果文件，计算分数
    while(<$FileHandler>)
    {
		my($line) = $_;
		chomp($line);

		my ( $ReadName, $Coverage ) = split( /\t/, $line, 2 );
        
        if ( exists( $ReadScoreHash{$ReadName} ) ) {
            $ReadScoreHash{$ReadName} -> [$PSLFileCount] = $Coverage;
        }
		
	} #while end

    $PSLFileCount ++;

}

# 关闭文件
foreach my $FileHandler ( @PSLFileArray ) {
  close $FileHandler;
}

open(OUTFILE, ">$FastAOutFile") || die "cannot open file $FastAOutFile \n" ;

foreach my $ReadName ( sort { $a cmp $b } keys %ReadScoreHash ) {

    print OUTFILE $ReadName."\t";
    print OUTFILE join( "\t", @{$ReadScoreHash{$ReadName}} );
    print OUTFILE "\n";

}

close(OUTFILE);
