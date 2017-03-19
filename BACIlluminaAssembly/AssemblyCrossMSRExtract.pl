use strict;
use warnings;
use Cwd;

# Maryland CA批量组装

my $Core  = 8;
my $Kmer  = 70;
my $Queue = 'vvl';

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

&createJob( );

sub createJob{

    my $CurrentDirectory = getcwd;
        
    if(! -e   "AssemblyCrossMSRResult"){
        mkdir "AssemblyCrossMSRResult";
    }
    
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

        foreach my $SampleVertical ( @SampleArrayVertical ){

            system("cp ./AssemblyCrossMSR/$SampleHorizontal-$SampleVertical/CA/10-gapclose/genome.ctg.fasta ./AssemblyCrossMSRResult/$SampleHorizontal-$SampleVertical.ctg.fasta");
        
        } # SampleVertical

	} # SampleHorizontal



} # createJob


