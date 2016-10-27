#!/usr/bin/perl
###########################################
# Todo: convert vcf to annovar acceptable 
#		annotation database.
# Usage: cat *.vcf|perl prepare_annovar.pl
# Author: Xiao-yu Zuo 
# Date: 2016-08-25
###########################################
use strict;
my $cohort="CANCER_FREE"; # specify the cohort name
my ($chr,$pos,$ref,$alt,$info);
# AC=#(mutant alleles); AF=AC/#(covered alleles)=AC/(2*NC) ; NC=#(mutant samples); NF=NC/#(covered samples); COVF=#(covered samples)/#(total samples)= 1 - sample_missing_rate
print "#Chr\tStart\tEnd\tRef\tAlt\tSYSUCC_${cohort}_AC\tSYSUCC_${cohort}_AF\tSYSUCC_${cohort}_NC\tSYSUCC_${cohort}_NF\tSYSUCC_${cohort}_COVF\n";
while(<>){
	chomp;
	next if /^#/;
	my @lines=split /\t/;
	($chr,$pos,$ref,$alt)=@lines[(0,1,3,4)];
	$chr=~s/^chr//;
	#print $chr,"\t",$pos,"\t",$ref,"\t",$alt,"\n";
	my @GENO;
	foreach my $gen (@lines[9..$#lines]){
		my @tmp=split /:/,$gen;
		push @GENO,$tmp[0];
	}
	my @alt_arr=split /,/,$alt;
	my $N=$#alt_arr+1;
	## if there is multiple ALT alleles in one loci (MNP or indel), split it and report seperately
	for(my $i=1;$i<=$N;$i++){
		next if $alt_arr[$i-1] eq "<*:DEL>"; # skip strange allele <*:DEL>
		my $ac=0; # allele count
		my $nc=0; # sample count
		my $nf=0; # allele freq
		my $af=0; # sample freq
		my $cc=0; # covered count
		my $cf=0; # covered freq
		foreach my $gen (@GENO){
			next if($gen eq "./."); # skip if missing (uncovered)
			$cc ++;
			my $t=grep{$_ eq $i} split /\//,$gen;
			$ac += $t;
			$nc ++ if $t > 0;
		}
		# if no covered samples, just set af, nf to 0
		if($cc > 0){
			$af=0.5*$ac/$cc;
			$nf=$nc/$cc; # only compute for covered samples
		}else{
			$af=0;
			$nf=0;
		}
		$cf=$cc/($#GENO+1);
		my ($start,$end,$ref,$alt)=recode_variant($pos,$ref,$alt_arr[$i-1]);
		print join "\t",($chr,$start,$end,$ref,$alt,$ac,$af,$nc,$nf,$cf);
		print "\n";
	}
}

### recode start, end, ref, alt based on its status: SNV, INSERTION, DELETION
sub recode_variant(){
	my ($pos,$ref,$alt)=@_;
	my ($start,$end);
	my $len_ref=length($ref);
	my $len_alt=length($alt);
	my $shift=0; # shift for pos, only for indel
	my ($b_ref,$b_alt); # current base for ref and alt
	if($len_ref==$len_alt){ #SNV
		if($len_alt==1){
			return ($pos,$pos,$ref,$alt);
		}else{ # e.g. AG->GG
			$shift=1;
			my $t1=substr $ref, $shift+1, $len_alt-$shift-1;
			my $t2=substr $alt, $shift+1, $len_alt-$shift-1;
			$t1 eq $t2 or die "$ref and $alt is strange in $pos\n";
			$ref=substr $ref, 0, $shift;
			$alt=substr $alt, 0, $shift;
			$start=$pos;
			$end=$pos;
		}
	}elsif($len_ref > $len_alt){ #deletion
		# compute shift.
		for(my $i=0;$i<$len_alt;$i++){ 
			$b_ref=substr $ref, $i, 1;
			$b_alt=substr $alt, $i, 1;
			last if $b_ref ne $b_alt;
			$shift++ if $b_ref eq $b_alt;
		}
		$start=$pos+$shift;
		$end=$pos+$len_ref-1;
		$alt=$shift == $len_alt? "-" : substr $alt, $shift;
		$ref=substr $ref, $shift;
		# if($shift == $len_alt){
		# 	$alt="-";
		# 	$ref=substr $ref, $shift;
		# }
		# if($shift<$len_alt){
		# 	my $t1=substr $ref, $shift+1, $len_alt-$shift-1;
		# 	my $t2=substr $alt, $shift+1, $len_alt-$shift-1;
		# 	$t1 eq $t2 or die "$ref and $alt is strange in $pos\n";
		# 	$alt="-";
		# 	$ref=sb
		# }
	}else{ #insertion
		# compute shift.
		for(my $i=0;$i<$len_ref;$i++){ 
			$b_ref=substr $ref, $i, 1;
			$b_alt=substr $alt, $i, 1;
			last if $b_ref ne $b_alt;
			$shift++ if $b_ref eq $b_alt;
		}
		$start=$pos+$shift;
		$end=$start;
		$ref=$shift == $len_ref? "-" : substr $ref, $shift;
		$alt=substr $alt, $shift;		
	}
	return ($start,$end,$ref,$alt);
}
