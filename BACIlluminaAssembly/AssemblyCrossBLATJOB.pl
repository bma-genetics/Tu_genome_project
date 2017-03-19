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

if(! -e 'AssemblyCrossBLAT'){
	mkdir 'AssemblyCrossBLAT';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
	
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {
        foreach my $SampleVertical ( @SampleArrayVertical ){
        
            open OUT,">$SampleHorizontal-$SampleVertical-AssemblyCrossBLAT".".pbs";
            print OUT <<SET;
#PBS -N CrossBLAT$SampleHorizontal$SampleVertical
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/

time blat ./AssemblyCopy/$SampleHorizontal.fasta \\
		  ./FilteredFasta/$SampleVertical-Filtered.fasta \\
		  ./AssemblyCrossBLAT/$SampleVertical-$SampleHorizontal.psl \\
		  -minIdentity=98 \\
		  -minScore=90 \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -out=psl \\
		  -maxGap=0            

SET
            close OUT;

        # 用来控制提交的JOBS数
        my $Jobs=`qstat -u bma | wc -l`;

        while( $Jobs > 1050 )
        {
            print "JOB Remain $Jobs\n";
            sleep 30;
            print `date`;
            $Jobs=`qstat -u bma | wc -l`;
        }

		system("echo $SampleHorizontal-$SampleVertical-AssemblyCrossBLAT".".pbs");
        system("qsub $SampleHorizontal-$SampleVertical-AssemblyCrossBLAT".".pbs");
        
        } # Vertical
        
	} # Horizontal

#==============================================================================================

    foreach my $SampleVertical ( @SampleArrayVertical ) {
        foreach my $SampleHorizontal ( @SampleArrayHorizontal ){
        
            open OUT,">$SampleVertical-$SampleHorizontal-AssemblyCrossBLAT".".pbs";
            print OUT <<SET;
#PBS -N CrossBLAT$SampleVertical$SampleHorizontal
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/

time blat ./AssemblyCopy/$SampleVertical.fasta \\
		  ./FilteredFasta/$SampleHorizontal-Filtered.fasta \\
		  ./AssemblyCrossBLAT/$SampleHorizontal-$SampleVertical.psl \\
		  -minIdentity=98 \\
		  -minScore=90 \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -out=psl \\
		  -maxGap=0           

SET
            close OUT;

        # 用来控制提交的JOBS数
        my $Jobs=`qstat -u bma | wc -l`;

        while( $Jobs > 1050 )
        {
            print "JOB Remain $Jobs\n";
            sleep 30;
            print `date`;
            $Jobs=`qstat -u bma | wc -l`;
        }

		system("echo $SampleVertical-$SampleHorizontal-AssemblyCrossBLAT".".pbs");
        system("qsub $SampleVertical-$SampleHorizontal-AssemblyCrossBLAT".".pbs");
        
        } # Vertical
	} # Horizontal
    

}
