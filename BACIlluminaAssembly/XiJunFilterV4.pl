use strict;
use warnings;
use Bio::SeqIO;

# 找出单端在载体上的序列V4，不切Reads
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

# 读入文件1的Read名字
open(INFILE,  "$PSLFileR1") || die "cannot open file $PSLFileR1 in DirectCount\n" ;

    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
		
		next if( $line !~ /^(\d+)\t(\d+)\t(\d+)/);
		
		my @Coulumn = split(/\t/, $line);
		
		# Match + Mismatch 大于90% 长度的Reads扔掉
		if ( ( $Coulumn[0] + $Coulumn[1] )/$Coulumn[10] >= $Threshold ) {
			$NameHashR1{$Coulumn[9]} = 0;
		}
		
	} #while end

close(INFILE) ;

# 从文件2找到单端在载体上的Read名字
open(INFILE,  "$PSLFileR2") || die "cannot open file $PSLFileR2 in DirectCount\n" ;

    while(<INFILE>)
    {
		my($line) = $_;
		chomp($line);
	
		next if( $line !~ /^(\d+)\t(\d+)\t(\d+)/);
		
		my @Coulumn = split(/\t/, $line);
		
		# Match + Mismatch 大于90% 长度的Reads扔掉
		if ( ( $Coulumn[0] + $Coulumn[1] )/$Coulumn[10] >= $Threshold ) {
			$NameHashR2{$Coulumn[9]} = 0;
		}
		
	} #while end

close(INFILE) ;

# 过滤掉单端在载体上的Reads

# R1
my $InFileR1  = Bio::SeqIO->new(-file => "$FastAInFileR1" , '-format' => 'Fasta');
my $OutFileR1 = Bio::SeqIO->new(-file => ">$FastAOutFileR1" , '-format' => 'Fasta');

while ( my $seq = $InFileR1 -> next_seq() ) {

	my $Scaffold = $seq->id;
	
	# 被记录的Reads
	if ( exists ( $NameHashR1{ $Scaffold } ) ) {

		# 如果是90%以上载体的Reads直接扔掉
		next if ( $NameHashR1{ $Scaffold } == 0 );	

	# 没被记录的Read直接输出
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
	
	# 被记录的Reads
	if ( exists ( $NameHashR2{ $Scaffold } ) ) {
		
		# 如果是90%以上载体的Reads直接扔掉
		next if ( $NameHashR2{ $Scaffold } == 0 );	
		
	# 没被记录的Read直接输出
	}else {

		$OutFileR2 -> width(500);
		$OutFileR2 -> write_seq($seq);
		
	}

}