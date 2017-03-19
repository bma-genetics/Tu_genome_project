use strict;
use warnings;
use Bio::SeqIO;

# ����һ������е�Reads�ڽ������г��ֵĸ��Ƕ�

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

# �ж�����������������ڴ�ֱ������ˮƽ����
if ( grep {/$PreFixQuery/} @SampleArrayVertical ) {
    
    # ��ֱ���򣬴�ˮƽ�����ļ�
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {
        # PSL�ȶ��ļ�
        my ($PSLFile) = "./AssemblyCrossCount/$PreFixQuery-$SampleHorizontal.txt";
        # print $PSLFile."\n";
        open my $FileHandler, "<$PSLFile";
        push ( @PSLFileArray, $FileHandler );

    }
    
}elsif ( grep {/$PreFixQuery/} @SampleArrayHorizontal ) {
    
    # ˮƽ���򣬴򿪴�ֱ�����ļ�
    foreach my $SampleVertical ( @SampleArrayVertical ) {
        # PSL�ȶ��ļ�
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

# ��ȡReadName,������Ӧ����������
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
    
    # ��ȡCount����ļ����������
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

# �ر��ļ�
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
