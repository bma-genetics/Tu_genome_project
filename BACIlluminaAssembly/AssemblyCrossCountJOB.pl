use strict;
use warnings;
use Cwd;

# 计算一个混池中的Reads出现次数

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

my @ArrayIndex = ( \@SampleArrayHorizontal, \@SampleArrayVertical );

if(! -e 'AssemblyCrossCount'){
	mkdir 'AssemblyCrossCount';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
	
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

            open OUT,">$SampleHorizontal-AssemblyCrossCount".".pbs";
            print OUT <<SET;
#PBS -N CrossCount$SampleHorizontal
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $Queue

cd \$PBS_O_WORKDIR

SET
            close OUT;
            
        foreach my $SampleVertical ( @SampleArrayVertical ){
            
            open OUT,">>$SampleHorizontal-AssemblyCrossCount".".pbs";
            print OUT <<SET;

time perl AssemblyCrossCount.pl \\
          ./AssemblyCrossBLAT/$SampleHorizontal-$SampleVertical.psl \\
          ./AssemblyCrossCount/$SampleHorizontal-$SampleVertical.txt

time perl AssemblyCrossCount.pl \\
          ./AssemblyCrossBLAT/$SampleVertical-$SampleHorizontal.psl \\
          ./AssemblyCrossCount/$SampleVertical-$SampleHorizontal.txt

SET
            close OUT;
            
        } # SampleVertical
        
        system("echo $SampleHorizontal-AssemblyCrossCount".".pbs");
        system("qsub $SampleHorizontal-AssemblyCrossCount".".pbs");

    } # SampleHorizontal

}

