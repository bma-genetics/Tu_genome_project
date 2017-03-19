use strict;
use warnings;
use Bio::SeqIO;

# ����һ������е�Reads�Խ���صĸ��Ƕ�

my $InFile  = shift;
my $OutFile = shift;

my %ReadsNameHash = ();

open( INFILE, "<$InFile") || die "cannot open file $InFile \n" ;

    # ��ȡPSL�ȶԽ���ļ����������
    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
        
		next if( $line !~ /^(\d+)\t(\d+)\t(\d+)/);
		
		my @Coulumn = split( /\t/, $line, 14 );
		
        my $Match    = $Coulumn[0];
        my $MisMatch = $Coulumn[1];
        my $GAPS     = $Coulumn[2].$Coulumn[3].$Coulumn[4].$Coulumn[5].$Coulumn[6].$Coulumn[7];
        my $QGAP     = $Coulumn[5];
        my $TGAP     = $Coulumn[7];
        my $ReadName = $Coulumn[9];
        my $QSize    = $Coulumn[10];
        my $QStart   = $Coulumn[11];
        my $QEnd     = $Coulumn[12];
        
        # ȡ��ÿ��Reads�����ĸ��Ƕ�
        if ( ! exists( $ReadsNameHash{$ReadName} ) ) {
            $ReadsNameHash{$ReadName} = sprintf ( "%.2f", ( $Match / $QSize ) );
        }elsif ( ( $Match / $QSize ) > $ReadsNameHash{$ReadName} ) {
            $ReadsNameHash{$ReadName} = sprintf ( "%.2f", ( $Match / $QSize ) );
		}
		
	} #while end

close( INFILE );

open( OUTFILE, ">$OutFile") || die "cannot open file $OutFile \n" ;
    
    foreach my $ReadName ( sort { $a cmp $b } keys %ReadsNameHash ) {
        print OUTFILE $ReadName."\t";
        print OUTFILE $ReadsNameHash{$ReadName}."\n";
    }

close( OUTFILE );