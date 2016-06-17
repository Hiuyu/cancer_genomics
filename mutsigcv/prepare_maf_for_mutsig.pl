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
		if($effect eq "unknown"){
			print STDERR "\"unknown\" effect found in <\"",join "\t", @arr,"\">, this annotation will be skipped.\n";
#print STDERR $_,"\n";
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
		$effect="noncoding" if grep{$_ =~/\b$func\b/} @noncoding_arr;
		$effect="silent" if grep{$_ =~/\b$func\b/} @silent_arr;
		$effect="nonsilent" if grep{$_ =~ /\b$func\b/} @nonsilent_arr;
		$effect="null" if grep{$_ =~ /$func/} @null_arr;
	}
	return $effect;	
}


sub check_and_split_annot{
	my $annot=shift;
	my $index=shift;
	my @split_annot;
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
			# do nothing
		}else{
			if($gene_func eq "intergenic"){ 
				$gene=deal_with_intergenic($gene, $arr_annot[$header{'GeneDetail_refGene'}]);
			}else{
				#warn "$gene\n";
				$gene=deal_with_multigene($gene);
				#die "$gene\n";
			}
#warn "$_\n";
#warn "-----------------------$gene--------------------------\n";
		}
	}else{
		($gene,$gene_func,$exonic_func)=deal_with_multifunc($gene,$gene_func,$exonic_func);
warn "$_\n";
warn "-----------------------$gene --- $gene_func --- $exonic_func --------------------------\n";
	}	
	$arr_annot[$index->[0]]=$gene;
	$arr_annot[$index->[1]]=$gene_func;
	$arr_annot[$index->[2]]=$exonic_func;
	push @split_annot, \@arr_annot;
	return @split_annot;
}
	
### pick the nearest gene, if Func_refGene == "intergenic"
sub deal_with_intergenic{
	my ($ref_gene, $dist_str)=@_;
	my $gene;	
	my $min_pos;
	if($dist_str=~/dist=(\d+);dist=(\d+)/){
		if($1 eq "NONE" && $2 eq "NONE"){
			warn "WARNING: $ref_gene, $dist_str has no dist value\n";
			return ".";
		}elsif($1 eq "NONE"){
			$min_pos=1;
		}elsif($2 eq "NONE"){
			$min_pos=0;
		}else{
			$min_pos=$1 < $2 ? 0 : 1;
		}
	}
	my @tmp=split /,/, $ref_gene;
	$gene=$tmp[$min_pos];
	return $gene;
}

### get the first non "LOC" gene
sub deal_with_multigene{
	my @genes=split /,/, shift;
	my $gene="";
	foreach my $g (@genes){
		$gene=$g;
		if($gene=~/^LOC/){
			next;
		}else{
			last;
		}
	}
	return $gene;
}

### get the most severous annotation, if has multiple Func_refGene
sub deal_with_multifunc{
	# gene, func_refGene, Func_exonic
	my @fg=split /;/, $_[0];
 	my @g=split /,/,$_[1];
        my @fe=split /;/,$_[2];
	my %fg_piority=(
		"exonic" => 1,
		"splicing" => 2,
		"UTR3" => 3,
		"UTR5" => 3,
		"intronic" => 4,
		"ncRNA_exonic" => 6,
		"ncRNA_intronic" => 7,
		"upstream" => 5,
		"downstream" => 5,
		"intergenic" => 8
		);
	my @fg_weight=map{$fg_piority{$_}} split /;/, $_[0];
	my ($min,$index)=index_min(@fg_weight);
	if($min==1){
		return ($g[$index], $fg[$index], $_[2]);
	}else{
		return ($g[$index], $fg[$index], ".");
	}
}

### find the min and index from an array
sub index_min{
	my ($min,$index)=(-1,-1);
	for(my $i=0;$i<@_;$i++){
		($min,$index)=($_[$i],0) if $index<0;
		if($min>$_[$i]){
			$min=$_[$i];
			$index=$i;
		}
	}
	return ($min, $index);
}
