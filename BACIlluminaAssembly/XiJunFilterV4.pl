use strict;
use warnings;
use Bio::SeqIO;

# �ҳ������������ϵ�����V4������Reads
my $PreFixR1  = shift;
my $PreFixR2  = shift;

my $Threshold = shift;

my $PSLFileR1  = './XiJunPSLData/'.$PreFixR1.'.psl';
my $PSLFileR2  = './XiJunPSLData/'.$PreFixR2.'.psl';

my $FastAInFileR1  = './BarCodeFasta/'.$PreFixR1.'.fasta';
my $FastAInFileR2  = './BarCodeFasta/'.$PreFixR2.'.fasta';

my $FastAOutFileR1  = './XiJunFilteredFasta/'.$PreFixR1.'-XiJunFiltered.fasta';
my $FastAOutFileR2  = './XiJunFilteredFasta/'.$PreFixR2.'-XiJunFiltered.fasta';

my %NameHashR1 = ();
my %NameHashR2 = ();

# �����ļ�1��Read����
open(INFILE,  "$PSLFileR1") || die "cannot open file $PSLFileR1 in DirectCount\n" ;

    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
		
		next if( $line !~ /^(\d+)\t(\d+)\t(\d+)/);
		
		my @Coulumn = split(/\t/, $line);
		
		# Match + Mismatch ����90% ���ȵ�Reads�ӵ�
		if ( ( $Coulumn[0] + $Coulumn[1] )/$Coulumn[10] >= $Threshold ) {
			$NameHashR1{$Coulumn[9]} = 0;
		}
		
	} #while end

close(INFILE) ;

# ���ļ�2�ҵ������������ϵ�Read����
open(INFILE,  "$PSLFileR2") || die "cannot open file $PSLFileR2 in DirectCount\n" ;

    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
	
		next if( $line !~ /^(\d+)\t(\d+)\t(\d+)/);
		
		my @Coulumn = split(/\t/, $line);
		
		# Match + Mismatch ����90% ���ȵ�Reads�ӵ�
		if ( ( $Coulumn[0] + $Coulumn[1] )/$Coulumn[10] >= $Threshold ) {
			$NameHashR2{$Coulumn[9]} = 0;
		}
		
	} #while end

close(INFILE) ;

# ���˵������������ϵ�Reads

# R1
my $InFileR1  = Bio::SeqIO->new(-file => "$FastAInFileR1" , '-format' => 'Fasta');
my $OutFileR1 = Bio::SeqIO->new(-file => ">$FastAOutFileR1" , '-format' => 'Fasta');

while ( my $seq = $InFileR1 -> next_seq() ) {

	my $Scaffold = $seq->id;
	
	# ����¼��Reads
	if ( exists ( $NameHashR1{ $Scaffold } ) ) {

		# �����90%���������Readsֱ���ӵ�
		next if ( $NameHashR1{ $Scaffold } == 0 );	

	# û����¼��Readֱ�����
	}else {
	
		$OutFileR1 -> width(500);
		$OutFileR1 -> write_seq($seq);
	
	}

}

# R2
my $InFileR2  = Bio::SeqIO->new(-file => "$FastAInFileR2" , '-format' => 'Fasta');
my $OutFileR2 = Bio::SeqIO->new(-file => ">$FastAOutFileR2" , '-format' => 'Fasta');

while ( my $seq = $InFileR2 -> next_seq() ) {

	my $Scaffold = $seq->id;
	
	# ����¼��Reads
	if ( exists ( $NameHashR2{ $Scaffold } ) ) {
		
		# �����90%���������Readsֱ���ӵ�
		next if ( $NameHashR2{ $Scaffold } == 0 );	
		
	# û����¼��Readֱ�����
	}else {

		$OutFileR2 -> width(500);
		$OutFileR2 -> write_seq($seq);
		
	}

}