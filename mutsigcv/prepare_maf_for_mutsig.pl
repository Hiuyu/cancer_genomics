#!/usr/bin/perl
#############################
# Title: prepare mutsig input file from annovar output
# Author: Xiao-yu Zuo
# Date: 2016-06-16
#############################
use strict;

my %header;
my $sample=shift @ARGV;
my $INPUT="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-annovar_output/${sample}.hg19_multianno.txt";
my @require_header=qw/Gene_refGene Chr Start Ref Alt/;
my @sub_arr;
open FILE, "<", $INPUT or die "$INPUT not found\n";
while(<FILE>){
	chomp;
	if($.==1){
		%header=get_header($_);
		print "gene\tchr\tstart\tref_allele\tnewbase\teffect\tpatient\n";
		next;
	}
	my @arr=split /\t/;
	foreach(@require_header){
		push @sub_arr, $arr[$header{$_}];
	}
	my $effect = define_effect($arr[$header{'Func_refGene'}], $arr[$header{'ExonicFunc_refGene'}]);
	push @sub_arr, $effect;
	push @sub_arr, $sample;
	








	undef @sub_arr;
}

sub get_header{
	my %name;
	my $line=$_;
	my @arr=split /\t/, $_; 
	my $count=0;
	foreach my $item (@arr){
		$name{$item}=$count;
		$count++;
	}
	return %name;
}	


sub define_effect{
	my $func_gene=shift;
	my $func_exonic=shift;
	my $func;
	my $effect="unknown";
	my @noncoding_arr=qw/intergenic intronic ncRNA_exonic ncRNA_intronic upstream downstream UTR3 UTR5 splicing/;
	my @silent_arr=qw/synonymous_SNV/;
	my @nonsilent_arr=qw/nonsynonymous_SNV stopgain stoplost/;
	my @null_arr=qw/nonframeshift_deletion nonframeshift_insertion frameshift_insertion frameshift_deletion/;
	
	if ($func_exonic=~/unknown/){
		return "unknown";
	}else{
		$func_exonic=~s/\s+/_/g;
		$func=($func_exonic eq '.' && $func_gene ne 'exonic') ? $func_gene : $func_exonic;
		$effect="noncoding" if grep{$_ eq $func} @noncoding_arr;
		$effect="silent" if grep{$_ eq $func} @silent_arr;
		$effect="nonsilent" if grep{$_ eq $func} @nonsilent_arr;
		$effect="null" if grep{$_ eq $func} @null_arr;
	}
	return $effect;	
}
