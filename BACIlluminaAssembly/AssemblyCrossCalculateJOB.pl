
use strict;
use warnings;
use Cwd;

my $Core  = 1;
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

if(! -e 'AssemblyCrossCalculate'){
    mkdir 'AssemblyCrossCalculate';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
    
    foreach my $VerticalSample ( @SampleArrayVertical ) {

            open OUT,">$VerticalSample-AssemblyCrossCalculate".".pbs";
            print OUT <<SET;
#PBS -N AssemblyCrossCalculate$VerticalSample
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR

date

SET
            close OUT;  
            
        foreach my $HorizontalSample ( @SampleArrayHorizontal ) {
            
            open OUT,">>$VerticalSample-AssemblyCrossCalculate".".pbs";
            print OUT <<SET;


time python coveragev_easy_pslv3.py \\
            -target ./AssemblyCopy/$VerticalSample.fasta \\
            -query ./FilteredFasta/$HorizontalSample-Filtered.fasta \\
            -psl ./AssemblyCrossBLAT/$VerticalSample-$HorizontalSample.psl \\
            -out ./AssemblyCrossCalculate/out.txt \\
            -depth ./AssemblyCrossCalculate/depth.txt \\
            > ./AssemblyCrossCalculate/$VerticalSample-$HorizontalSample.txt

time python coveragev_easy_pslv3.py \\
            -target ./AssemblyCopy/$HorizontalSample.fasta \\
            -query ./FilteredFasta/$VerticalSample-Filtered.fasta \\
            -psl ./AssemblyCrossBLAT/$HorizontalSample-$VerticalSample.psl \\
            -out ./AssemblyCrossCalculate/out.txt \\
            -depth ./AssemblyCrossCalculate/depth.txt \\
            > ./AssemblyCrossCalculate/$HorizontalSample-$VerticalSample.txt

date

SET
            close OUT;

        } # SampleArrayHorizontal
        
        system("echo $VerticalSample-AssemblyCrossCalculate".".pbs");
        # system("qsub $VerticalSample-AssemblyCrossCalculate".".pbs");

    }# SampleArrayVertical


} # createJob

