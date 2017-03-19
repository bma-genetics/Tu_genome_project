use strict;
use warnings;
use Cwd;

# ��BLAT�ȶ��г��Reads�ŵ���Ӧ��BAC��ȥ

my $HorizontalFile = 'Horizontal.txt';
my $VerticalFile   = 'Vertical.txt';
my $Queue          = 'vvl';

my @SampleArrayHorizontal = ();
my @SampleArrayVertical   = ();

my $ExtractQuota    = 1;    # ���Cell��
my $MinimumCoverage = 0.90; # ��С����
my $Times           = 2;    # ����

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

if(! -e 'AssemblyCrossExtract'){
	mkdir 'AssemblyCrossExtract';
}

&createJob( );

sub createJob{
    
    my $CurrentDirectory = getcwd;
	
    foreach my $SampleArray ( @ArrayIndex ) {
        
        foreach my $Sample ( @{$SampleArray} ){
            
            open OUT,">$Sample-AssemblyCrossExtract".".pbs";
            print OUT <<SET;
#PBS -N CrossExtract$Sample
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q $Queue

cd \$PBS_O_WORKDIR

#    perl AssemblyCrossExtract.pl <Sample> <ExtractQuota> <MinimumCoverage> <Times>
time perl AssemblyCrossExtract.pl $Sample $ExtractQuota $MinimumCoverage $Times > $Sample.log

SET
            close OUT;

            system("echo $Sample-AssemblyCrossExtract".".pbs");
            system("qsub $Sample-AssemblyCrossExtract".".pbs");
        
        } # SampleArrayHorizontal

    } # SampleArrayVertical

}

