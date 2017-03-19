use strict;
use warnings;
use Cwd;

# MatePair连接BAC片段

my $HorizontalFile = 'Horizontal.txt';
my $VerticalFile   = 'Vertical.txt';
my $Queue          = 'high';
my $Core           = 1;

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

if(! -e 'AssemblyMatePair'){
	mkdir 'AssemblyMatePair';
}

if(! -e 'AssemblyMatePairResult'){
	mkdir 'AssemblyMatePairResult';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
	
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {
        
        foreach my $SampleVertical ( @SampleArrayVertical ){
            
            open OUT,">$SampleHorizontal-$SampleVertical-AssemblyMatePair".".pbs";
            print OUT <<SET;
#PBS -N MatePair$SampleHorizontal-$SampleVertical
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/AssemblyMatePair

hostname

time perl /public/share/bma/SSPACE/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl \\
          -l ../libraries.txt \\
          -s ../AssemblyCrossMSRResult/$SampleHorizontal-$SampleVertical.ctg.fasta \\
          -T $Core \\
          -k 5 \\
          -a 0.7 \\
          -x 0 \\
          -b $SampleHorizontal-$SampleVertical

time rm ./$SampleHorizontal-$SampleVertical/reads/*.fa
time rm ./$SampleHorizontal-$SampleVertical/alignoutput/*
time rm ./$SampleHorizontal-$SampleVertical/intermediate_results/*
time rm ./$SampleHorizontal-$SampleVertical/pairinfo/*

time cp ./$SampleHorizontal-$SampleVertical/$SampleHorizontal-$SampleVertical.final.scaffolds.fasta ../AssemblyMatePairResult/$SampleHorizontal-$SampleVertical.final.scaffolds.fasta
          
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
            
            system("echo $SampleHorizontal-$SampleVertical-AssemblyMatePair".".pbs");
            system("qsub $SampleHorizontal-$SampleVertical-AssemblyMatePair".".pbs");
        
        } # SampleArrayHorizontal

    } # SampleArrayVertical

}

