#!/usr/bin/perl -w
use strict;

my @tempArr;
my @tempArr1;
my $input;

##read Ta gff
$input = "./MCScanTuTa/ta.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %ta_locus = ();
my %ta_pos = ();
while(<INFILE>){    
   chomp;
   if($_ =~ /^\#/){next;}
   @tempArr = split(/\s+/,$_);
   $ta_locus{$tempArr[2]}{$tempArr[3]}= $tempArr[1];
   $ta_pos{$tempArr[1]}= $tempArr[2]."\t".$tempArr[3];
}
close INFILE;
#for my $key (sort {$a<=>$b} keys %ta_locus) {
#   print $key."\t".$ta_locus{$key}."\n";
#}
#exit;

##read gff Tu.
$input = "./MCScanTuTa/tu.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %tu_locus = ();
my %tu_pos = ();
while(<INFILE>){                   
   chomp;
   @tempArr = split(/\s+/,$_);
   if($tempArr[0] eq "tu3"){
       $tu_locus{$tempArr[2]}{$tempArr[3]}= $tempArr[1];
       $tu_pos{$tempArr[1]}= $tempArr[2]."\t".$tempArr[3];
   }
}
close INFILE;
#for my $key (sort {$a<=>$b} keys %tu_locus) {
#    for my $skey (sort {$a<=>$b} keys %{$tu_locus{$key}}){
#	print $key."\t".$skey."\t".$tu_locus{$key}{$skey}."\n";
#    }
#}
#exit;

##read gff bd.
$input = "./MCScanTaBd/bd.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   

my %bd_pos = ();
while(<INFILE>){                   
   chomp;
   @tempArr = split(/\s+/,$_);
   if($tempArr[0] eq "bd2"){
       $bd_pos{$tempArr[1]}= $tempArr[2]."\t".$tempArr[3];
   }
}
close INFILE;

#for my $key (sort keys %bd_pos) {
#   print $key."\t".$bd_pos{$key}."\n";
#}
#exit;

##read gff rice
$input = "./MCScanTaOs/os.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %os_pos = ();
while(<INFILE>){                   
   chomp;
   @tempArr = split(/\s+/,$_);
   if($tempArr[0] eq "os1"){
       $os_pos{$tempArr[1]}= $tempArr[2]."\t".$tempArr[3];
   }
}
close INFILE;
#for my $key (sort keys %os_pos) {
#   print $key."\t".$os_pos{$key}."\n";
#}
#exit;
##read Ta Tu collinearity

##read gff sb
$input = "./MCScanTaSb/tasb.gff";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %sb_pos = ();
while(<INFILE>){                   
   chomp;
   @tempArr = split(/\s+/,$_);
   if($tempArr[0] eq "sb3"){
       $sb_pos{$tempArr[1]}= $tempArr[2]."\t".$tempArr[3];
   }
}
close INFILE;

##read no the bases of collinear blocks
$input = "./MCScanTuTa/tuta.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %tatu_col = ();
my %tuta_col = ();
my %block_dic = ();
my $block_id = 0;
my $tu_start = "NA";
my $tu_end ;
my $ta_start ;
my $ta_end;

while(<INFILE>){
   if($_ =~ /^\#/){
       if($_ =~ /Alignment/){
	   @tempArr = split(/\s+/,$_);
	   if($tu_start ne "NA"){
	       #print $tu_start."\n";
	       @tempArr1 = split(/\s+/,$tu_pos{$tu_start});
	       my @tempArr2 = split(/\s+/,$tu_pos{$tu_end});
	       if($tempArr1[0] < $tempArr2[0]){
		   $block_dic{$block_id}= $tu_start."\t".$tu_end;
	       }
	       else{
		   $block_dic{$block_id}= $tu_end."\t".$tu_start;
	       }
	       @tempArr1 = split(/\s+/,$ta_pos{$ta_start});
	       @tempArr2 = split(/\s+/,$ta_pos{$ta_end});
	       if($tempArr1[0] < $tempArr2[0]){
		   $block_dic{$block_id} .= "\t".$ta_start."\t".$ta_end;
	       }
	       else{
		   $block_dic{$block_id} .= "\t".$ta_end."\t".$ta_start;
	       }
	   }
	   $block_id = substr($tempArr[2],0,-1);
	   $tu_start = "NA";
       }
       next;
   }
   if($_ =~ /TuG1812G03/){
       $_ =~ s/^\s+//;                
       chomp;
       @tempArr = split(/\s+/,$_);
       $tempArr[0] =~ s/\-//;
       $tempArr[1] =~ s/\://;
       #$tatu_col{$tempArr[3]} = $tempArr[4]."\t".$tempArr[1]."\t".$tempArr[2];
       #$tuta_col{$tempArr[4]} = $tempArr[3]."\t".$tempArr[1]."\t".$tempArr[2];
       $tatu_col{$block_id}{$tempArr[2]} = $tempArr[3];
       $tuta_col{$block_id}{$tempArr[3]} = $tempArr[2];
       #print $block_id."\n";
       if($tu_start eq "NA"){
	   $tu_start = $tempArr[3];
	   $ta_start = $tempArr[2];
       }
       else{
	   $tu_end = $tempArr[3];
	   $ta_end = $tempArr[2];
       }
   }
}
close INFILE;

#for my $key (sort {$a<=>$b} keys %block_dic) {
#   print $key."\t".$block_dic{$key}."\n";
#}
#exit;


##read bd Ta collinearity
$input = "./MCScanTaBd/tabd.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %bdta_col = ();
my %tabd_col = ();

while(<INFILE>){     
   if($_ =~ /^\#/){next;} 
   if($_ =~ /Bradi2g/){
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_);
       #$tempArr[0] =~ s/\-//;
       #$tempArr[1] =~ s/\://;
       if($tempArr[2] =~ /^Bradi/){
	   $bdta_col{$tempArr[2]} = $tempArr[3];
	   $tabd_col{$tempArr[3]} = $tempArr[2];
       }
       elsif($tempArr[2] =~ /^TRAE/){
	   print $_."\n";exit;
	   $bdta_col{$tempArr[1]} = $tempArr[2];
	   $tabd_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;
#for my $key (keys %bdta_col) {
#    print $key."\t".$bdta_col{$key}."\n";
#}
#exit;

##read bd Tu collinearity
$input = "./MCScanTuBd/tubd.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %bdtu_col = ();
my %tubd_col = ();

while(<INFILE>){     
   if($_ =~ /^\#/){next;}
   if($_ =~ /Bradi2g/ && $_ =~ /TuG1812G01/){
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_);
       $tempArr[0] =~ s/\-//;
       $tempArr[1] =~ s/\://;
       
       if($tempArr[2] =~ /^Bradi/){
	   $bdtu_col{$tempArr[2]} = $tempArr[3];
	   $tubd_col{$tempArr[3]} = $tempArr[2];
       }
       else{
	   $bdtu_col{$tempArr[1]} = $tempArr[2];
	   $tubd_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;
#for my $key (keys %bdtu_col) {
#    print $key."\t".$bdtu_col{$key}."\n";
#}
#exit;

##read Tu Os collinearity
$input = "./MCScanTuOs/tuos.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %ostu_col = ();
my %tuos_col = ();

while(<INFILE>){  
   if($_ =~ /^\#/){next;} 
   if($_ =~ /Os01g/){
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_);
       $tempArr[0] =~ s/\-//;
       $tempArr[1] =~ s/\://;
       if($tempArr[2] =~ /^Os01g/){
	   $ostu_col{$tempArr[2]} = $tempArr[3];
	   $tuos_col{$tempArr[3]} = $tempArr[2];
       }
       else{
	   $ostu_col{$tempArr[1]} = $tempArr[2];
	   $tuos_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;

#for my $key (keys %ostu_col) {
#    print $key."\t".$ostu_col{$key}."\n";
#}
#exit;

##read Ta Os collinearity
$input = "./MCScanTaOs/taos.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %osta_col = ();
my %taos_col = ();

while(<INFILE>){
   if($_ =~ /^\#/){next;}  
   if($_ =~ /Os01g/){                 
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_); 
       if($tempArr[2] =~ /^Os01g/){
	   $osta_col{$tempArr[2]} = $tempArr[3];
	   $taos_col{$tempArr[3]} = $tempArr[2];
       }
       else{
	   $osta_col{$tempArr[1]} = $tempArr[2];
	   $taos_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;

##read Tu sb collinearity
$input = "./MCScanTuSb/tusb.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %sbtu_col = ();
my %tusb_col = ();

while(<INFILE>){  
   if($_ =~ /^\#/){next;} 
   if($_ =~ /Sobic.003G/){
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_);
       $tempArr[0] =~ s/\-//;
       $tempArr[1] =~ s/\://;
       if($tempArr[2] =~ /^Sobic.003G/){
	   $sbtu_col{$tempArr[2]} = $tempArr[3];
	   $tusb_col{$tempArr[3]} = $tempArr[2];
       }
       else{
	   $sbtu_col{$tempArr[1]} = $tempArr[2];
	   $tusb_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;

#for my $key (keys %sbtu_col) {
#    print $key."\t".$sbtu_col{$key}."\n";
#}
#exit;

##read Ta sb collinearity
$input = "./MCScanTaSb/tasb.collinearity";
open(INFILE, $input) ||  
   die "Unable to open $input for reading: $!\n";
   
my %sbta_col = ();
my %tasb_col = ();

while(<INFILE>){
   if($_ =~ /^\#/){next;}  
   if($_ =~ /Sobic.003G/){                 
       chomp;
       $_ =~ s/^\s+//; 
       @tempArr = split(/\s+/,$_); 
       if($tempArr[2] =~ /^Sobic.003G/){
	   $sbta_col{$tempArr[2]} = $tempArr[3];
	   $tasb_col{$tempArr[3]} = $tempArr[2];
       }
       else{
	   $sbta_col{$tempArr[1]} = $tempArr[2];
	   $tasb_col{$tempArr[2]} = $tempArr[1];
       }
   }
}
close INFILE;
#for my $key (keys %sbta_col) {
#    print $key."\t".$sbta_col{$key}."\n";
#}
#exit;
##Using Tu as reference
my $output = "Tu_ref_collinearity.out";
open(OUTFILE, ">".$output) ||
    die "Unable to open $output for writing: $!\n";

print OUTFILE "tuID\ttustart\ttuend\ttaID\ttastart\ttaend\tbdID\tbdstart\tbdend\tosID\tosstart\tosend\tsbID\tsbstart\tsbend\n";
for my $bkey (sort {$a<=>$b} keys %block_dic ){
    print OUTFILE "## Alignment ".$bkey.":\n";
    @tempArr = split(/\t/,$block_dic{$bkey});
    @tempArr1 = split(/\t/,$tu_pos{$tempArr[0]});
    $tu_start = $tempArr1[0];
    @tempArr1 = split(/\t/,$tu_pos{$tempArr[1]});
    $tu_end = $tempArr1[0];
    for my $key ( sort {$a<=>$b} keys %tu_locus) {
	if($key < $tu_start){next;}
	if($key > $tu_end){last;}
	for my $skey ( sort {$a<=>$b} keys %{$tu_locus{$key}}) {
	    my $tuid = $tu_locus{$key}{$skey};
	    print OUTFILE $tuid."\t".$tu_pos{$tuid};
	    if(exists $tuta_col{$bkey}{$tuid}){
		print OUTFILE "\t".$tuta_col{$bkey}{$tuid}."\t".$ta_pos{$tuta_col{$bkey}{$tuid}};
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check bd
	    if(exists $tubd_col{$tuid}){
		print OUTFILE "\t".$tubd_col{$tuid}."\t".$bd_pos{$tubd_col{$tuid}};
	    }
	    elsif(exists $tuta_col{$bkey}{$tuid}){
		if(exists $tabd_col{$tuta_col{$bkey}{$tuid}}){
		    print OUTFILE "\t".$tabd_col{$tuta_col{$bkey}{$tuid}};
		    if(!exists $bd_pos{$tabd_col{$tuta_col{$bkey}{$tuid}}}){
			print $tuid."\t".$tuta_col{$bkey}{$tuid}."\t".$tabd_col{$tuta_col{$bkey}{$tuid}}."\n";
			exit;
		    }
		    print OUTFILE "\t".$bd_pos{$tabd_col{$tuta_col{$bkey}{$tuid}}};
		}
		else{
		    print OUTFILE "\t--\t--\t--";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check os
	    if(exists $tuos_col{$tuid}){
		print OUTFILE "\t".$tuos_col{$tuid}."\t".$os_pos{$tuos_col{$tuid}};
	    }
	    elsif(exists $tuta_col{$bkey}{$tuid}){
		if(exists $taos_col{$tuta_col{$bkey}{$tuid}}){
		    print OUTFILE "\t".$taos_col{$tuta_col{$bkey}{$tuid}}."\t".$os_pos{$taos_col{$tuta_col{$bkey}{$tuid}}};
		}
		else{
		    print OUTFILE "\t--\t--\t--";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check sb
	    if(exists $tusb_col{$tuid}){
		print OUTFILE "\t".$tusb_col{$tuid}."\t".$sb_pos{$tusb_col{$tuid}}."\n";
	    }
	    elsif(exists $tuta_col{$bkey}{$tuid}){
		if(exists $tasb_col{$tuta_col{$bkey}{$tuid}}){
		    print OUTFILE "\t".$tasb_col{$tuta_col{$bkey}{$tuid}}."\t".$sb_pos{$tasb_col{$tuta_col{$bkey}{$tuid}}}."\n";
	    }
		else{
		    print OUTFILE "\t--\t--\t--\n";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--\n";
	    }
	}
    }
}
close OUTFILE;


##Using Ta as reference
$output = "Ta_ref_collinearity.out";
open(OUTFILE, ">".$output) ||
    die "Unable to open $output for writing: $!\n";
print OUTFILE "taID\ttastart\ttaend\ttuID\ttustart\ttuend\tbdID\tbdstart\tbdend\tosID\tosstart\tosend\tsbID\tsbstart\tsbend\n";

for my $bkey (sort {$a<=>$b} keys %block_dic ){
    print OUTFILE "## Alignment ".$bkey.":\n";
    @tempArr = split(/\t/,$block_dic{$bkey});
    @tempArr1 = split(/\t/,$ta_pos{$tempArr[2]});
    $ta_start = $tempArr1[0];
    @tempArr1 = split(/\t/,$ta_pos{$tempArr[3]});
    $ta_end = $tempArr1[0];
    for my $key ( sort {$a<=>$b} keys %ta_locus) {
	if($key < $ta_start){next;}
	if($key > $ta_end){last;}
	for my $skey ( sort {$a<=>$b} keys %{$ta_locus{$key}}) {
	    my $taid = $ta_locus{$key}{$skey};
	    print OUTFILE $taid."\t".$ta_pos{$taid};
	    if(exists $tatu_col{$bkey}{$taid}){
		if(!exists $tu_pos{$tatu_col{$bkey}{$taid}}){
		    print $tatu_col{$bkey}{$taid}."\n";
		}
		print OUTFILE "\t".$tatu_col{$bkey}{$taid}."\t".$tu_pos{$tatu_col{$bkey}{$taid}};
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check bd
	    if(exists $tabd_col{$taid}){
		print OUTFILE "\t".$tabd_col{$taid}."\t".$bd_pos{$tabd_col{$taid}};
	    }
	    elsif(exists $tatu_col{$bkey}{$taid}){
		if(exists $tubd_col{$tatu_col{$bkey}{$taid}}){
		    print OUTFILE "\t".$tubd_col{$tatu_col{$bkey}{$taid}}."\t".$bd_pos{$tubd_col{$tatu_col{$bkey}{$taid}}};
		}
		else{
		    print OUTFILE "\t--\t--\t--";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check os
	    if(exists $taos_col{$taid}){
		print OUTFILE "\t".$taos_col{$taid}."\t".$os_pos{$taos_col{$taid}};
	    }
	    elsif(exists $tuta_col{$bkey}{$taid}){
		if(exists $tuos_col{$tuta_col{$bkey}{$taid}}){
		    print OUTFILE "\t".$tuos_col{$tatu_col{$bkey}{$taid}}."\t".$os_pos{$tuos_col{$tatu_col{$bkey}{$taid}}};
		}
		else{
		    print OUTFILE "\t--\t--\t--";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--";
	    }
	    ##check sb
	    if(exists $tasb_col{$taid}){
		print OUTFILE "\t".$tasb_col{$taid}."\t".$sb_pos{$tasb_col{$taid}}."\n";
	    }
	    elsif(exists $tuta_col{$bkey}{$taid}){
		if(exists $tusb_col{$tuta_col{$bkey}{$taid}}){
		    print OUTFILE "\t".$tusb_col{$tatu_col{$bkey}{$taid}}."\t".$sb_pos{$tusb_col{$tatu_col{$bkey}{$taid}}}."\n";
		}
		else{
		    print OUTFILE "\t--\t--\t--\n";
		}
	    }
	    else{
		print OUTFILE "\t--\t--\t--\n";
	    }
	}
    }
}
close OUTFILE;

##Using bd as reference
#$output = "Bd_ref_collinearity.out";
#open(OUTFILE, ">".$output) ||
#    die "Unable to open $output for writing: $!\n";
#print "bdID\tbdstart\tbdend\taID\ttastart\ttaend\ttuID\ttustart\ttuend\n";
#for my $key ( sort {$a<=>$b} keys %bd_locus) {
#    @tempArr = split(/\t/,$bd_locus{$key});
#    print OUTFILE $taid."\t".$key."\t".$tempArr[1];
#    if(exists $bdta_col{$taid}){
#	print OUTFILE "\t".$bdta_col{$taid}."\t".$ta_pos{$bdta_col{$taid}};
#    }
#    else{
#	print OUTFILE "\t--\t--\t--";
#    }
#    if(exists $bdtu_col{$taid}){
#	print OUTFILE "\t".$bdtu_col{$taid}."\t".$tu_pos{$bdtu_col{$taid}}."\n";
#    }
#    else{
#	print OUTFILE "\t--\t--\t--\n";
#    }
#}
#close OUTFILE;

print "end of format_collinear.pl\n";
