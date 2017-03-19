use strict;
use warnings;
use Cwd;
use Bio::SeqIO;

# 把抽出来的Reads分成双端数据，交叉的BAC混在一起组装

my $SampleHorizontal    = shift;
my $SampleVertical      = shift;
my $FastqPathHorizontal = shift;
my $FastqPathVertical   = shift;
my $SampleH             = shift;
my $SampleV             = shift;
my $Query               = '';
my $LineCount           = 0;
my $OutputFlag          = 0;

# 用来记录Reads是否取出
my %ReadPairHash   = ();

# 从两个方向上的cell里取出ReadName
open(INFILE,  "./AssemblyCrossExtract/$SampleVertical-$SampleHorizontal.txt") || die "cannot open file \n" ;

    while(<INFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        my ( $ReadName, $Others ) = split ( /\t/, $line, 2 );
        $ReadName =~ s/\-1// if ( $ReadName =~ /\-1/ );
        $ReadName =~ s/\-2// if ( $ReadName =~ /\-2/ );
		$ReadPairHash{$ReadName} = 0;
        
    } #while end

close(INFILE) ;

open(INFILE,  "./AssemblyCrossExtract/$SampleHorizontal-$SampleVertical.txt") || die "cannot open file \n" ;

    while(<INFILE>)
    {
        my($line) = $_;
        chomp($line);
        
        my ( $ReadName, $Others ) = split ( /\t/, $line, 2 );
        $ReadName =~ s/\-1// if ( $ReadName =~ /\-1/ );
        $ReadName =~ s/\-2// if ( $ReadName =~ /\-2/ );
		$ReadPairHash{$ReadName} = 0;
        
    } #while end

close(INFILE) ;

open( OUTFILER1,  ">./AssemblyCrossReads/$SampleHorizontal-$SampleVertical-R1.fastq") || die "cannot open file \n" ;
open( OUTFILER2,  ">./AssemblyCrossReads/$SampleHorizontal-$SampleVertical-R2.fastq") || die "cannot open file \n" ;

# ==========================================================================================================

$Query = $SampleHorizontal;
$Query =~ s/$SampleH/Query/;

# open( INSOURCEFILER1,  "<$FastqPathHorizontal/$Query/R1.fastq") || die "cannot open file \n" ;
open( INSOURCEFILER1, "gunzip -c $FastqPathHorizontal/$Query/R1.fastq.gz |") || die "can't open pipe to $FastqPathHorizontal/$Query/R1.fastq.gz" ;

    $LineCount  = 0;
    $OutputFlag = 0;
    
    while(<INSOURCEFILER1>)
    {
		my($line) = $_;
		chomp($line);
		
		$LineCount++;
		
		# 取Fastq文件的第一行名称
		if( $LineCount % 4 == 1 ) {
			
			if ( $line =~ m/^@(\S+)\s\S+/ ) {
				
				my $ReadName = $1;
				
				if ( exists( $ReadPairHash{$ReadName} ) ) {
					print OUTFILER1 $line."\n";
					$OutputFlag = 1;
				}
				
			}
			
		} # end of if

		# 输出后几行数据
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 2 && $OutputFlag == 1 );
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 3 && $OutputFlag == 1 );
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 0 && $OutputFlag == 1 );
		# 清空FLAG
		$OutputFlag = 0 if( $LineCount % 4 == 0 && $OutputFlag == 1 );
	}

close(INSOURCEFILER1) ;

# open( INSOURCEFILER2,  "<$FastqPathHorizontal/$Query/R2.fastq") || die "cannot open file \n" ;
open( INSOURCEFILER2, "gunzip -c $FastqPathHorizontal/$Query/R2.fastq.gz |") || die "can't open pipe to $FastqPathHorizontal/$Query/R2.fastq.gz" ;

    $LineCount  = 0;
    $OutputFlag = 0;
    
    while(<INSOURCEFILER2>)
    {
		my($line) = $_;
		chomp($line);
		
		$LineCount++;
		
		# 取Fastq文件的第一行名称
		if( $LineCount % 4 == 1 ) {
			
			if ( $line =~ m/^@(\S+)\s\S+/ ) {
				
				my $ReadName = $1;
				
				if ( exists( $ReadPairHash{$ReadName} ) ) {
					print OUTFILER2 $line."\n";
					$OutputFlag = 1;
				}
				
			}
			
		} # end of if

		# 输出后几行数据
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 2 && $OutputFlag == 1 );
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 3 && $OutputFlag == 1 );
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 0 && $OutputFlag == 1 );
		# 清空FLAG
		$OutputFlag = 0 if( $LineCount % 4 == 0 && $OutputFlag == 1 );
	}

close(INSOURCEFILER2) ;

# ==========================================================================================================

$Query = $SampleVertical;
$Query =~ s/$SampleV/Query/;

# open( INSOURCEFILER1,  "<$FastqPathVertical/$Query/R1.fastq") || die "cannot open file \n" ;
open( INSOURCEFILER1, "gunzip -c $FastqPathVertical/$Query/R1.fastq.gz |") || die "can't open pipe to $FastqPathVertical/$Query/R1.fastq.gz" ;

    $LineCount  = 0;
    $OutputFlag = 0;
    
    while(<INSOURCEFILER1>)
    {
		my($line) = $_;
		chomp($line);
		
		$LineCount++;
		
		# 取Fastq文件的第一行名称
		if( $LineCount % 4 == 1 ) {
			
			if ( $line =~ m/^@(\S+)\s\S+/ ) {
				
				my $ReadName = $1;
				
				if ( exists( $ReadPairHash{$ReadName} ) ) {
					print OUTFILER1 $line."\n";
					$OutputFlag = 1;
				}
				
			}
			
		} # end of if

		# 输出后几行数据
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 2 && $OutputFlag == 1 );
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 3 && $OutputFlag == 1 );
		print OUTFILER1 $line."\n" if( $LineCount % 4 == 0 && $OutputFlag == 1 );
		# 清空FLAG
		$OutputFlag = 0 if( $LineCount % 4 == 0 && $OutputFlag == 1 );
	}

close(INSOURCEFILER1) ;

# open( INSOURCEFILER2,  "<$FastqPathVertical/$Query/R2.fastq") || die "cannot open file \n" ;
open( INSOURCEFILER2, "gunzip -c $FastqPathVertical/$Query/R2.fastq.gz |") || die "can't open pipe to $FastqPathVertical/$Query/R2.fastq.gz" ;

    $LineCount  = 0;
    $OutputFlag = 0;
    
    while(<INSOURCEFILER2>)
    {
		my($line) = $_;
		chomp($line);
		
		$LineCount++;
		
		# 取Fastq文件的第一行名称
		if( $LineCount % 4 == 1 ) {
			
			if ( $line =~ m/^@(\S+)\s\S+/ ) {
				
				my $ReadName = $1;
				
				if ( exists( $ReadPairHash{$ReadName} ) ){
					print OUTFILER2 $line."\n";
					$OutputFlag = 1;
				}
				
			}
			
		} # end of if

		# 输出后几行数据
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 2 && $OutputFlag == 1 );
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 3 && $OutputFlag == 1 );
		print OUTFILER2 $line."\n" if( $LineCount % 4 == 0 && $OutputFlag == 1 );
		# 清空FLAG
		$OutputFlag = 0 if( $LineCount % 4 == 0 && $OutputFlag == 1 );
	}

close(INSOURCEFILER2) ;

close(OUTFILER1) ;
close(OUTFILER2) ;

# my $AllCount = 0;
# my $PairCount = 0;
# my $R1SingleCount = 0;
# my $R2SingleCount = 0;

# print "$SampleHorizontal-$SampleVertical\n";
# print "$SampleVertical-$SampleHorizontal\n";
# print "\$AllCount\t$AllCount\n";
# print "\$PairCount\t$PairCount\n";
# print "\$R1SingleCount\t$R1SingleCount\n";
# print "\$R2SingleCount\t$R2SingleCount\n";