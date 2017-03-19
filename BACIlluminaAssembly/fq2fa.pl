#!/usr/bin/perl
use strict;
my $i;

# open(INPUT_FILE,"<$ARGV[0]");
open(OUTPUT_FILE,">$ARGV[0]");

while(<STDIN>){
	if($i%4==0){
		my $ReadName = substr($_,1);
		$ReadName =~ s/ /-/;
		$ReadName =~ s/:N:.*//;
		print OUTPUT_FILE ">".$ReadName;
	}
	else{
		if($i%4==1){
			print OUTPUT_FILE $_;
		}
	}
	$i++;
}

# close(INPUT_FILE);
close(OUTPUT_FILE);
