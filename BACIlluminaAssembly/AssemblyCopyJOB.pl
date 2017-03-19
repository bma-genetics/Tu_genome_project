use strict;
use warnings;
use Cwd;

my $Core  = 1;
my $Queue = 'low';

my $AssemblyResultPathS = "SED-AssemblyResultPathS";
my $AssemblyResultPathL = "SED-AssemblyResultPathL";

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

if(! -e 'AssemblyCopy'){
	mkdir 'AssemblyCopy';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
    
    foreach my $Sample ( @SampleArrayHorizontal ) {
            
            my $SeqNum = $Sample;
            $SeqNum    =~ s/E1X//g if ( $SeqNum =~ m/E1X/ ) ;
            $SeqNum    =~ s/a6X//g if ( $SeqNum =~ m/a6X/ ) ;
            $SeqNum    =~ s/G62-X//g if ( $SeqNum =~ m/G62-X/ ) ;
            $SeqNum    =~ s/G63-X//g if ( $SeqNum =~ m/G63-X/ ) ;
            $SeqNum    =~ s/PV50-P//g if ( $SeqNum =~ m/PV50-P/ ) ;
            $SeqNum    =~ s/PV51-P//g if ( $SeqNum =~ m/PV51-P/ ) ;
            $SeqNum    =~ s/C1X//g if ( $SeqNum =~ m/C1X/ ) ;
            $SeqNum    =~ s/C2X//g if ( $SeqNum =~ m/C2X/ ) ;
            $SeqNum    =~ s/C3X//g if ( $SeqNum =~ m/C3X/ ) ;
            $SeqNum    =~ s/C4X//g if ( $SeqNum =~ m/C4X/ ) ;
            
            open OUT,">AssemblyCopy$Sample".".pbs";
            print OUT <<SET;
#PBS -N AssemblyCopy$Sample
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/AssemblyCopy

touch ./$Sample.fasta
cp $AssemblyResultPathL/09-Query$SeqNum.utg.fasta ./$Sample.fasta

SET
            close OUT;
		system("echo AssemblyCopy$Sample".".pbs");
        system("qsub AssemblyCopy$Sample".".pbs");
    }

    foreach my $Sample ( @SampleArrayVertical ) {
            
            my $SeqNum = $Sample;
            $SeqNum    =~ s/E1Y//g if ( $SeqNum =~ m/E1Y/ ) ;
            $SeqNum    =~ s/a6Y//g if ( $SeqNum =~ m/a6Y/ ) ;
            $SeqNum    =~ s/G62-Y//g if ( $SeqNum =~ m/G62-Y/ ) ;
            $SeqNum    =~ s/G63-Y//g if ( $SeqNum =~ m/G63-Y/ ) ;
            $SeqNum    =~ s/PV50-V//g if ( $SeqNum =~ m/PV50-V/ ) ;
            $SeqNum    =~ s/PV51-V//g if ( $SeqNum =~ m/PV51-V/ ) ;
            $SeqNum    =~ s/C1Y//g if ( $SeqNum =~ m/C1Y/ ) ;
            $SeqNum    =~ s/C2Y//g if ( $SeqNum =~ m/C2Y/ ) ;
            $SeqNum    =~ s/C3Y//g if ( $SeqNum =~ m/C3Y/ ) ;
            $SeqNum    =~ s/C4Y//g if ( $SeqNum =~ m/C4Y/ ) ;
            
            open OUT,">AssemblyCopy$Sample".".pbs";
            print OUT <<SET;
#PBS -N AssemblyCopy$Sample
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/AssemblyCopy

touch ./$Sample.fasta
cp $AssemblyResultPathL/09-Query$SeqNum.utg.fasta ./$Sample.fasta

SET
            close OUT;
		system("echo AssemblyCopy$Sample".".pbs");
        system("qsub AssemblyCopy$Sample".".pbs");
    }

}
