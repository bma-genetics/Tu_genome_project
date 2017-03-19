use strict;
use warnings;
use Cwd;

# 取出每个BAC的总长矩阵

my $BACFile       = "ctg-2.txt";
my %BACLengthHash = ();

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

        my @Columns = split( /\s+/, $line );
        my $BACName = $Columns[-1];
        my $Sum     = $Columns[-2];
        $BACLengthHash{$BACName} = $Sum;
        
    } # end of while

close( BACFILE );


foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

    foreach my $SampleVertical ( @SampleArrayVertical ){
        
        if ( exists( $BACLengthHash{"$SampleHorizontal-$SampleVertical"} ) ) {
            # print "$SampleHorizontal-$SampleVertical\t";
            print $BACLengthHash{"$SampleHorizontal-$SampleVertical"}."\t";
        }else {
            # print "$SampleHorizontal-$SampleVertical\t";
            print "0\t";        
        }
        
    } # SampleVertical
    print "\n";
} # SampleHorizontal
