use strict;
use warnings;
use Cwd;

# 取出每个BAC的总长矩阵

my $BACFile       = "ctg-2.txt";
my $BACReadFile   = "BACReadsCount.txt";
my %BACLengthHash = ();
my %BACReadHash   = ();

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

open( BACFILE, "<$BACFile") || die "cannot open file $BACFile \n" ;

    while(<BACFILE>)
    {
        my($line) = $_;
        chomp($line);

        # 13      13      1       436     12316   12316   72973   72973   72973   122261  10X10-10Y10
        my @Columns = split( /\s+/, $line );
        my $BACName = $Columns[-1];
        my $Info    = $line;
        $BACLengthHash{$BACName} = join( "\t" , @Columns );
        
    } # end of while

close( BACFILE );

open( BACREADFILE, "<$BACReadFile") || die "cannot open file $BACReadFile \n" ;

    while(<BACREADFILE>)
    {
        my($line) = $_;
        chomp($line);

        # 10X10-10Y10	8777439
        my @Columns = split( /\s+/, $line );
        my $BACName = $Columns[0];
        my $Reads   = $Columns[1];
        $BACReadHash{$BACName} = $Reads;
        
    } # end of while

close( BACREADFILE );

foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

    foreach my $SampleVertical ( @SampleArrayVertical ){
        
        if ( exists( $BACLengthHash{"$SampleHorizontal-$SampleVertical"} ) ) {
            print $BACLengthHash{"$SampleHorizontal-$SampleVertical"}."\t";
        }else {
            print "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t$SampleHorizontal-$SampleVertical\t";
        }
      
        if ( exists( $BACReadHash{"$SampleHorizontal-$SampleVertical"} ) ) {
            print $BACReadHash{"$SampleHorizontal-$SampleVertical"}."\n";
        }else {
            print "0\n";
        }

    } # SampleVertical

} # SampleHorizontal
