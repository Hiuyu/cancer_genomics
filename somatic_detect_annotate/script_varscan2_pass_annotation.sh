#!/bin/bash
##################################
# Title:combine individual varscan2 output,
#       make annotation, and make a similar table like that from mutect
# Author: Xiao-yu Zuo
# Date: 2016-09-23
##################################
# grep ^# 3-mutect2_PASS/WES00003_pass_mutect2.vcf|perl -lane 'next if /^##GATK/ || /^##SAMPLE/; if(/^#CHROM/){print "##INFO=<ID=SAMPLE,Number=.,Type=String,Description=\"Sample ID\">";} print $_; ' > 3.1-mutect2_PASS/header.vcf
#########################
# Step 1: combine varscan2 snv output to a big vcf
#         remove chrM and other small contigs
# Assume the header file is copied here.
#########################
workdir="$PWD" # working path, all required files should be manually copied here
source_single_path="$PWD/../7-varscan_PASS" # the path of individual pass mutect2/mutect file.
OUT_PREFIX="varscan2_pass" # prefix for the final output
snv_surfix="varscan2_PASS_SNV.txt" # the surfix of individual pass vcf, used to get the ${sample}_$single_surfix, e.g. WES00014_mutect_PASS.vcf
indel_surfix="varscan2_PASS_INDEL.txt" # the surfix of individual pass vcf, used to get the ${sample}_$single_surfix, e.g. WES00014_mutect_PASS.vcf
PASS_INPUT="varscan2_pass.avinput"
sample_list="all_list.txt"


### Combine all file into a big one
cd $workdir
echo "combining varscan outputs"
test -f body.${OUT_PREFIX}.out && rm body.${OUT_PREFIX}.out
# if not append SAMPLE INFO, just add it 
for sample in `cat $sample_list`
do 
  echo $sample
  # formatting SNV output
  cat $source_single_path/${sample}_${snv_surfix} | perl -lwane '
    next unless $.>1 || $F[0] =~ /^chr[MXY\d][\d]*$/;
    $gt0="$F[2]/$F[2]";
    $gt1=$F[8]==0?"$F[2]/$F[3]":"$F[3]/$F[3]";
    # the first 5 columns were repeated, in order to keep original 0-based coordination to make consistent to the GATK output
    print join "\t",($F[0],$F[1],$F[1],$F[2],$F[3],$F[0],$F[1],$F[1],$F[2],$F[3],"'$sample'",$gt0,"$F[4],$F[5]",$gt1,"$F[8],$F[9]");
  ' >> body.${OUT_PREFIX}.out
  # formatting INDEL output
  cat $source_single_path/${sample}_${indel_surfix} | perl -lwane '
    ## remove small contig
    next unless $.>1 || $F[0] =~ /^chr[MXY\d][\d]*$/; 
    $chr = $F[0];
    $flag = substr($F[3],0,1); # get the type flag "-" => DEL, "+" => INS
    $seq = substr($F[3],1); # get the sequence of INS or DEL
    $up1 = $F[2]; # get the 0-based former allele
    ## define start, end
    $start1 = $F[1] + 1; # 1-based start position
    $start0 = $F[1]; # 0-based start position
    $end1 = $flag eq "-" ? $start1 + length($seq) -1 : $start1; # DEL => start + len(del_seq), INS => start
    $end0 = $flag eq "-" ? $start0 + length($seq)  : $start0; 
    ## define ref, alt
    $ref1 = "";
    $alt1 = "";
    $ref0 = ""; # ref allele in 0-based cooridnation
    $alt0 = "";
    $gt0N = ""; # genotype of NORMAL in 1-based coordinaiton
    $gt0T = "";
    if($flag eq "-"){ 
      # DEL
      $ref1 = $seq;
      $ref0 = "$up1$seq";
      $alt1 = "-";
      $alt0 = "$up1";
    }elsif($flag eq "+"){
      # INS
      $ref1 = "-";
      $ref0 = $up1;
      $alt1 = $seq;
      $alt0 = "$up1$seq";
    }else{
      die "$. have unknown flag: < $flag > in  $F[3]\n";
    }
    $gt0N = "$ref0/$ref0"; # use 0-based genotype
    $gt0T = $F[8] == 0 ? "$alt0/$alt0" : "$ref0/$alt0" ; 
    # INDEL should keep original 0-based coordination to make consistent to the GATK output
    print join "\t",($chr,$start1,$end1,$ref1,$alt1,$chr,$start0,$end0,$ref0,$alt0,"'$sample'",$gt0N,"$F[4],$F[5]",$gt0T,"$F[8],$F[9]");
  ' >> body.${OUT_PREFIX}.out
done
# sort and remove chrM and/or other contigs(if have?)
cat body.${OUT_PREFIX}.out|sort -V -k 1,1 -k 2,2n > $PASS_INPUT
rm body.${OUT_PREFIX}.out

##########################
# Step 2: make annotation on the big vcf
# 
##########################
ANV_PREFIX=$OUT_PREFIX
ANV_TABLE="$ANV_PREFIX.tsv"
fasta="$HOME/resource/DNASeq/hg19/ucsc_hg19.fa"
nt=8

echo "making annotation"
perl /data/home/zuoxy/apps/annovar/2015Dec14/table_annovar.pl $PASS_INPUT /data/home/zuoxy/resource/annovar_hg19/humandb/ -buildver hg19 -out $ANV_PREFIX -remove -protocol refGene,avsnp144,cosmic74,dbnsfp30a,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_sas,exac03,esp6500siv2_all,dbnsfp31a_interpro,dbscsnv11,kaviar_20150923,hrcr1,sysucc_cancer_free_841sample_20161027 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . --dot2underline -thread $nt

##########################
# Step 3: merge table
# 
##########################
echo "to table"
echo -e "CHROM\tStart\tEnd\tREF\tALT\tSAMPLE\tNORMAL.GT\tNORMAL.AD\tTUMOR.GT\tTUMOR.AD" > $ANV_PREFIX.part1.txt
cat $PASS_INPUT | cut -f 6-  >> $ANV_PREFIX.part1.txt
cat $ANV_PREFIX.hg19_multianno.txt|cut -f 6- > $ANV_PREFIX.part2.txt
paste $ANV_PREFIX.part1.txt $ANV_PREFIX.part2.txt  > ${ANV_PREFIX}_annotate.tsv
rm $ANV_PREFIX.part1.txt $ANV_PREFIX.part2.txt


