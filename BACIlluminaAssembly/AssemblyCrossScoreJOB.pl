use strict;
use warnings;
use Cwd;

# 计算一个混池中的Reads在交叉混池中出现的次数

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

if(! -e 'AssemblyCrossScore'){
	mkdir 'AssemblyCrossScore';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
	
    foreach my $SampleArray ( @ArrayIndex ) {
        
        foreach my $Sample ( @{$SampleArray} ){
            
            open OUT,">$Sample-AssemblyCrossScore".".pbs";
            print OUT <<SET;
#PBS -N CrossScore$Sample
#PBS -j oe
#PBS -l nodes=1:ppn=8
#PBS -q $Queue

cd \$PBS_O_WORKDIR

time perl AssemblyCrossScore.pl $Sample

SET
            close OUT;

            system("echo $Sample-AssemblyCrossScore".".pbs");
            system("qsub $Sample-AssemblyCrossScore".".pbs");
        
        } # Sample

    } # SampleArray

}

