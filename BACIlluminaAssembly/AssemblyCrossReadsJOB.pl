use strict;
use warnings;
use Cwd;

my $FastqPathHorizontal = "SED-QueryFastqPathS";
my $FastqPathVertical   = "SED-QueryFastqPathL";
my $SampleH             = 'SED-SampleS';
my $SampleV             = 'SED-SampleL';

my $HorizontalFile = 'Horizontal.txt';
my $VerticalFile   = 'Vertical.txt';
my $Queue          = 'vvl';

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

if(! -e 'AssemblyCrossReads'){
	mkdir 'AssemblyCrossReads';
}

&createJob( );

sub createJob{

    my $CurrentDirectory = getcwd;
	
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

            open OUT,">$SampleHorizontal-AssemblyCrossReads".".pbs";
            print OUT <<SET;
#PBS -N CrossReads$SampleHorizontal
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $Queue

cd \$PBS_O_WORKDIR/

SET
            close OUT;
			
        foreach my $SampleVertical ( @SampleArrayVertical ){
            
            open OUT,">>$SampleHorizontal-AssemblyCrossReads".".pbs";
            print OUT <<SET;

time perl AssemblyCrossReads.pl $SampleHorizontal $SampleVertical $FastqPathHorizontal $FastqPathVertical $SampleH $SampleV
		 
SET
            close OUT;
           
        } # Vertical
		
		system("echo $SampleHorizontal-AssemblyCrossReads".".pbs");
        system("qsub $SampleHorizontal-AssemblyCrossReads".".pbs");

	} # Horizontal


}
