use warnings ;
use strict ;
use Bio::SeqIO;

my $InputFile  = shift; #'TUG.E0470.J07-bac.fasta';
my $OutputFile = shift; #'TUG.E0470.J07.scf.fasta';
#my $min_len=shift;

my $in  = Bio::SeqIO->new(-file => "$InputFile" ,    '-format' => 'Fasta');
my $out = Bio::SeqIO->new(-file => ">$OutputFile" ,  '-format' => 'Fasta');

my %LengthHash = ();
my %SeqHash    = ();

while ( my $seq = $in->next_seq() ) {

    my $ScfID  = $seq -> id;
    my $Length = $seq -> length;
    my $Seq    = $seq -> seq;

    $LengthHash{$ScfID} = $Length;
    $SeqHash{$ScfID}    = $Seq;

}

my $ScaffoldNum = 1;
foreach my $ScfID ( sort { $LengthHash{$b} <=> $LengthHash{$a} } keys %LengthHash ) {
    # scaffold30|size478
    my $Length   = $LengthHash{$ScfID};
#    next if($Length<$min_len);
    my $SubSeq   = Bio::Seq->new( -seq => $SeqHash{$ScfID} , -display_id => "scaffold$ScaffoldNum|size$Length" );

    $out->width(100);
    $out->write_seq($SubSeq);

    $ScaffoldNum ++

}
