use strict;
use warnings;

# 混池分离-去细菌和载体

my @PoolReadsPosition = ( "SED-QueryFastqPathS",
                          "SED-QueryFastqPathL",
                         );

my $ZaitiRef          = "/public/share/bma/Reference-BWA/ZaiTi_Rice_500bp.fasta";
my $XiJunRef          = "/public/share/bma/WheatAssembly/XiJunAlignment/EscherichiaColi.fasta";
my $minIdentity       = 95;
my $minScore          = 50;
my $Threshold         = 0.5;
my $SampleS           = 'SED-SampleS';
my $SampleL           = 'SED-SampleL';

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

my @ArrayIndex = ( \@SampleArrayHorizontal, \@SampleArrayVertical );

if(! -e 'BarCodeFasta'){
	mkdir 'BarCodeFasta';
}

if(! -e 'ZaiTiPSLData'){
	mkdir 'ZaiTiPSLData';
}

if(! -e 'XiJunPSLData'){
	mkdir 'XiJunPSLData';
}

if(! -e 'ZaiTiFilteredFasta'){
	mkdir 'ZaiTiFilteredFasta';
}

if(! -e 'XiJunFilteredFasta'){
	mkdir 'XiJunFilteredFasta';
}

if(! -e 'ValidatePSLData'){
	mkdir 'ValidatePSLData';
}

if(! -e 'FilteredFasta'){
	mkdir 'FilteredFasta';
}

&createJob( );

sub createJob{
    
    my $Num = -1;
    
    foreach my $SampleArray ( @ArrayIndex ) {
        
        $Num ++;
        
        foreach my $Sample ( @{$SampleArray} ){

            my $Query = '';
            
            if ( $Sample =~ /$SampleS/ ) {
                $Query = $Sample;
                $Query =~ s/$SampleS/Query/;
            }elsif ( $Sample =~ /$SampleL/ ) {
                $Query = $Sample;
                $Query =~ s/$SampleL/Query/;
            }
            
            open OUT,">$Sample-fqfa".".pbs";
            print OUT <<SET;
#PBS -N $Sample-fqfa
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -q vvl

cd \$PBS_O_WORKDIR

zcat $PoolReadsPosition[$Num]/$Query/R1.fastq | perl fq2fa.pl  ./BarCodeFasta/$Sample-R1.fasta
zcat $PoolReadsPosition[$Num]/$Query/R2.fastq | perl fq2fa.pl  ./BarCodeFasta/$Sample-R2.fasta
# time perl fq2fa.pl $PoolReadsPosition[$Num]/$Query/R1.fastq ./BarCodeFasta/$Sample-R1.fasta
# time perl fq2fa.pl $PoolReadsPosition[$Num]/$Query/R2.fastq ./BarCodeFasta/$Sample-R2.fasta

time blat $ZaitiRef \\
		  ./BarCodeFasta/$Sample-R1.fasta \\
		  ./ZaiTiPSLData/$Sample-R1.psl \\
		  -minIdentity=$minIdentity \\
		  -minScore=$minScore \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -maxGap=0

time blat $ZaitiRef \\
		  ./BarCodeFasta/$Sample-R2.fasta \\
		  ./ZaiTiPSLData/$Sample-R2.psl \\
		  -minIdentity=$minIdentity \\
		  -minScore=$minScore \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -maxGap=0

time blat $XiJunRef \\
		  ./BarCodeFasta/$Sample-R1.fasta \\
		  ./XiJunPSLData/$Sample-R1.psl \\
		  -minIdentity=$minIdentity \\
		  -minScore=$minScore \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -maxGap=0

time blat $XiJunRef \\
		  ./BarCodeFasta/$Sample-R2.fasta \\
		  ./XiJunPSLData/$Sample-R2.psl \\
		  -minIdentity=$minIdentity \\
		  -minScore=$minScore \\
		  -tileSize=18 \\
		  -stepSize=19 \\
          -fastMap \\
		  -maxGap=0

time perl XiJunFilterV4.pl $Sample-R1 $Sample-R2 $Threshold
time perl ZaiTiFilterV4.pl $Sample-R1 $Sample-R2 $Threshold

cat ./ZaiTiFilteredFasta/$Sample-R1-ZaiTiFiltered.fasta > ./FilteredFasta/$Sample-R1-Filtered.fasta
cat ./ZaiTiFilteredFasta/$Sample-R2-ZaiTiFiltered.fasta > ./FilteredFasta/$Sample-R2-Filtered.fasta
cat ./FilteredFasta/$Sample-R1-Filtered.fasta ./FilteredFasta/$Sample-R2-Filtered.fasta > ./FilteredFasta/$Sample-Filtered.fasta

# time blat $ZaitiRef \\
# 		  ./ZaiTiFilteredFasta/$Sample-R1-ZaiTiFiltered.fasta \\
# 		  ./ValidatePSLData/$Sample-R1-ZaiTiFiltered.psl \\
# 		  -minIdentity=$minIdentity \\
# 		  -minScore=$minScore \\
# 		  -tileSize=18 \\
# 		  -stepSize=19 \\
#         -fastMap \\
# 		  -maxGap=0
# 
# time blat $ZaitiRef \\
# 		  ./ZaiTiFilteredFasta/$Sample-R2-ZaiTiFiltered.fasta \\
# 		  ./ValidatePSLData/$Sample-R2-ZaiTiFiltered.psl \\
# 		  -minIdentity=$minIdentity \\
# 		  -minScore=$minScore \\
# 		  -tileSize=18 \\
# 		  -stepSize=19 \\
#         -fastMap \\
# 		  -maxGap=0
# 
#           
# time blat $XiJunRef \\
# 		  ./ZaiTiFilteredFasta/$Sample-R1-ZaiTiFiltered.fasta \\
# 		  ./ValidatePSLData/$Sample-R1-XiJunFiltered.psl \\
# 		  -minIdentity=$minIdentity \\
# 		  -minScore=$minScore \\
# 		  -tileSize=18 \\
# 		  -stepSize=19 \\
#         -fastMap \\
# 		  -maxGap=0
# 
# time blat $XiJunRef \\
# 		  ./ZaiTiFilteredFasta/$Sample-R2-ZaiTiFiltered.fasta \\
# 		  ./ValidatePSLData/$Sample-R2-XiJunFiltered.psl \\
# 		  -minIdentity=$minIdentity \\
# 		  -minScore=$minScore \\
# 		  -tileSize=18 \\
# 		  -stepSize=19 \\
#         -fastMap \\
# 		  -maxGap=0
  
SET
            close OUT;
            system("echo $Sample-fqfa".".pbs");
            system("qsub $Sample-fqfa".".pbs");
            
        }
}


}
