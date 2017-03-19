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
        
    if(! -e   "AssemblyCrossMSR"){
        mkdir "AssemblyCrossMSR";
    }
    
    foreach my $SampleHorizontal ( @SampleArrayHorizontal ) {

        foreach my $SampleVertical ( @SampleArrayVertical ){
            
            my $PathR1 = "$CurrentDirectory/AssemblyCrossReads/$SampleHorizontal-$SampleVertical-R1.fastq";
            my $PathR2 = "$CurrentDirectory/AssemblyCrossReads/$SampleHorizontal-$SampleVertical-R2.fastq";
            
            if(! -e   "./AssemblyCrossMSR/$SampleHorizontal-$SampleVertical"){
                mkdir "./AssemblyCrossMSR/$SampleHorizontal-$SampleVertical";
            }
            
            open OUT,">$CurrentDirectory/AssemblyCrossMSR/$SampleHorizontal-$SampleVertical/Wheat-maryland";
            print OUT <<SET;
#example configuration file for rhodobacter sphaeroides assembly from GAGE project (http://gage.cbcb.umd.edu)

#DATA is specified as type {PE,JUMP,OTHER}= two_letter_prefix mean stdev fastq(.gz)_fwd_reads fastq(.gz)_rev_reads
#NOTE that PE reads are always assumed to be  innies, i.e. --->  <---, and JUMP are assumed to be outties <---    --->; if there are any jump libraries that are innies, such as longjump, specify them as JUMP and specify NEGATIVE mean
#rev reads are optional for PE libraries and mandatory for JUMP libraries
#any OTHER sequence data (454, Sanger, Ion torrent, etc) must be first converted into Celera Assembler compatible .frg files (see http://wgs-assembler.sourceforge.com)

DATA

# InsertSize
PE= pe 280 42 $PathR1 $PathR2

END

PARAMETERS

#this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE=$Kmer

#set this to 1 for Illumina-only assemblies and to 0 if you have 2x or more long (Sanger, 454) reads
USE_LINKING_MATES=1

#this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and something large (300) for mammals
LIMIT_JUMP_COVERAGE = 100

#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. for mammals do not set cgwErrorRate above 0.15!!!
CA_PARAMETERS = merylMemory=8192 cgwErrorRate=0.15 ovlMemory=8GB

#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100
KMER_COUNT_THRESHOLD = 1

#auto-detected number of cpus to use
NUM_THREADS= $Core

#this is mandatory jellyfish hash size
JF_SIZE=500000000

#this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes
DO_HOMOPOLYMER_TRIM=0

END

SET
            close OUT;
            
            open OUT,">$CurrentDirectory/AssemblyCrossMSR/$SampleHorizontal-$SampleVertical/MSR$SampleHorizontal-$SampleVertical-Assembly.pbs";
            print OUT <<SET;
#PBS -N MSR_Wheat$SampleHorizontal-$SampleVertical
#PBS -j oe
#PBS -l nodes=1:ppn=$Core
#PBS -q $Queue

cd \$PBS_O_WORKDIR/

time assemble.sh
		 
SET
            close OUT;
		
        # 用来控制提交的JOBS数
        my $Jobs=`qstat -u bma | wc -l`;

        while( $Jobs > 500 )
        {
            print "JOB Remain $Jobs\n";
            sleep 30;
            print `date`;
            $Jobs=`qstat -u bma | wc -l`;
        }
        
		system("echo $CurrentDirectory/MSR$SampleHorizontal-$SampleVertical-Assembly.pbs");
        chdir "$CurrentDirectory/AssemblyCrossMSR/$SampleHorizontal-$SampleVertical";
        # system("/public/share/bma/MaSuRCA/MaSuRCA-2.2.1/bin/masurca Wheat-maryland");
        # sleep(5);
        system("qsub MSR$SampleHorizontal-$SampleVertical-Assembly.pbs");
        chdir "$CurrentDirectory";
        
        } # SampleVertical

	} # SampleHorizontal



} # createJob
