#!/usr/bin/perl
#####################################################################
# Title: convert format for IntOGen website(www.intogen.org/analysis)
#        from 'maf' to tab-separated format.
# Author: Xiao-yu Zuo
# Date: 2016-06-27
#####################################################################

### assume stdin
while(<>){
	next if $.==1;
	chomp;
	@F=split /\t/;
	$strand="+";
	$chr=$F[1];
	$start=$F[2];
	$end=$start;
	$ref=$F[3];
	$alt=$F[4];
	$sample=$F[16];
	## decide $start and $end ($start + length($ref) -1 )
	if($ref eq '-'){
	### insertion
	### start and end are reversed in the case of insertion
		$start=$end;
		$end=$end-1;
	}elsif($alt eq '-'){
	### deletion
	### let the $end plus the length of deletion
		$end=$start+length($ref)-1;
	}else{
	### SNV
	### do nothing
		1;
	}	
	$allele="$ref>$alt";
	## print output
	print "$chr\t$start\t$end\t$strand\t$allele\t$sample\n";
}
