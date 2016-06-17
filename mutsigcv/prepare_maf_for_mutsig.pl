#!/usr/bin/perl
#############################
# Title: prepare mutsig input file from annovar output
# Author: Xiao-yu Zuo
# Date: 2016-06-16
#############################
use strict;
use List::Util qw(min max reduce);

my %header;
my $sample=shift @ARGV;
my $INPUT="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-annovar_output/${sample}.hg19_multianno.txt";
my @require_header=qw/Gene_refGene Chr Start Ref Alt/;
my @sub_arr;
my @index_multianno_check;
open FILE, "<", $INPUT or die "$INPUT not found\n";
while(<FILE>){
	chomp;
	if($.==1){
		%header=get_header($_);
		@index_multianno_check=($header{'Gene_refGene'},$header{'Func_refGene'},$header{'ExonicFunc_refGene'});
		print "gene\tchr\tstart\tref_allele\tnewbase\teffect\tpatient\n";
		next;
	}
	### if multiple annotation, do it one by one
	my @arr_annot=check_and_split_annot($_, \@index_multianno_check);
	foreach my $line (@arr_annot){
		my @arr=@$line;
		foreach(@require_header){
			push @sub_arr, $arr[$header{$_}];
		}
		my $effect = define_effect($arr[$header{'Func_refGene'}], $arr[$header{'ExonicFunc_refGene'}]);
		### if found 'unknown' effect, warning and skip
		if($effect eq "null"){
			#print STDERR "\"unknown\" effect found in <\"",join "\t", @arr,"\">, this annotation will be skipped.\n";
print STDERR $_,"\n";
			undef @sub_arr;
			next;
		}
		push @sub_arr, $effect;
		push @sub_arr, $sample;
		print join "\t", @sub_arr,"\n";
		undef @sub_arr;
	}
		
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
	my $effect="null";
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


sub check_and_split_annot{
	my $annot=shift;
	my $index=shift;
	my @split_annot;
	my $is_split=0;
	my @arr_annot=split /\t/, $annot;
	my $gene=$arr_annot[$index->[0]];
	my $gene_func=$arr_annot[$index->[1]];
	my $exonic_func=$arr_annot[$index->[2]];
	### if found ",", split into multiple lines
	my $num_gene=split /,/,$gene;
	my $num_gene_func=split /;/,$gene_func;
	my $num_exonic_func=split /;/,$exonic_func;
	my @nums=($num_gene,$num_gene_func,$num_exonic_func);
	if($num_gene_func==1){
		if($num_gene==1 && $num_exonic_func==1){
			$is_split=0;
		}else{
			if($gene_func eq "intergenic"){ ######################## not finish ###########################
				
	### if func_gene==intergenic, do nothing
	if($arr_annot[$index->[1]] eq "intergenic"){
		$is_split=0;
	}
	if ($is_split){
		my $num_gene=split /,/,$arr_annot[$index->[0]];
		my $num_gene_func=split /,/,$arr_annot[$index->[1]];
		my $num_exonic_func=split /,/,$arr_annot[$index->[2]];
		
		my $i=0;
		while($i<$num_gene){
			my @temp;
			foreach my $j(0..@arr_annot){
				if(grep{$_ == $j} @$index){
					my @tmp_split=split /,/,$arr_annot[$j];
					push @temp,$tmp_split[$i];

				}else{
					push @temp,$arr_annot[$j];
				}
			}
			push @split_annot,\@temp;
			$i++;
		}
	}else{
		push @split_annot, \@arr_annot;
	}
	return @split_annot;
}
	


