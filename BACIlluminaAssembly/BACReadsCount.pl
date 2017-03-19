use warnings ;
use strict ;
use Bio::SeqIO;
use Bio::Perl;

my $Num = 1;

my $R1Results = `abyss-fac -t 0 ./AssemblyCrossReads/*X*-*Y*-R1.fastq | grep -v sum | awk '{print \$(NF-1)}'`;
my @R1ResultsArray = split( /\n/, $R1Results );

my $R2Results = `abyss-fac -t 0 ./AssemblyCrossReads/*X*-*Y*-R2.fastq | grep -v sum | awk '{print \$(NF-1)}'`;
my @R2ResultsArray = split( /\n/, $R2Results );

my $BACNameResults = `abyss-fac -t 0 ./AssemblyCrossReads/*X*-*Y*-R1.fastq | grep -v sum | awk '{print \$NF}' | awk -F/ '{print \$NF}' | sed 's/-R1.fastq//'`;
my @BACNameResultsArray = split( /\n/, $BACNameResults );

for ( my $ResultsIndex = 0; $ResultsIndex < scalar( @R1ResultsArray ); $ResultsIndex++ ) {
    print $BACNameResultsArray[$ResultsIndex]."\t";
    # print $R1ResultsArray[$ResultsIndex]."\t";
    # print $R2ResultsArray[$ResultsIndex]."\t";
    print $R1ResultsArray[$ResultsIndex]+$R2ResultsArray[$ResultsIndex];
    print "\n";
    $Num ++;
}

# $R1Results = `abyss-fac -t 0 ./AssemblyCrossReads/Query??/R1.fastq | grep -v sum | awk '{print \$9}'`;
# @R1ResultsArray = split( /\n/, $R1Results );

# $R2Results = `abyss-fac -t 0 ./AssemblyCrossReads/Query??/R2.fastq | grep -v sum | awk '{print \$9}'`;
# @R2ResultsArray = split( /\n/, $R2Results );

# for ( my $ResultsIndex = 0; $ResultsIndex < scalar( @R1ResultsArray ); $ResultsIndex++ ) {
    # print $Num."\t";
    # print $R1ResultsArray[$ResultsIndex]."\t";
    # print $R2ResultsArray[$ResultsIndex]."\t";
    # print $R1ResultsArray[$ResultsIndex]+$R2ResultsArray[$ResultsIndex];
    # print "\n";
    # $Num ++;
# }

# for ( my $QueryNum = 1; $QueryNum <= 48; $QueryNum ++ ) {
    
    # my $TotalLength = 0;
    
    # my $in  = Bio::SeqIO->new(-file => "./FilteredFastq/Query$QueryNum/R1.fastq" , '-format' => 'Fastq');

    # while ( my $seq = $in->next_seq() ) {

        # my $Length = $seq->length;
        
        # $TotalLength = $TotalLength + $Length;

    # }

    # $in  = Bio::SeqIO->new(-file => "./FilteredFastq/Query$QueryNum/R2.fastq" , '-format' => 'Fastq');

    # while ( my $seq = $in->next_seq() ) {

        # my $Length = $seq->length;
        
        # $TotalLength = $TotalLength + $Length;

    # }

    # print "Query$QueryNum\t$TotalLength\n"

# }