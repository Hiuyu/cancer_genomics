#!/bin/sh
module load java/1.8.0_60
##############################################
# Title: MutSigCV analysis based on multiple 
# 	single vcfs from Mutect
# Author: Xiao-yu Zuo
# Date: 2016-06-05 
##############################################

VCF_PATH="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene"
IDS_LIST=$1 # please input the full path of the IDS file.
##############################################
# convert mutect vcf output to maf (mutation 
# annotation format)
##############################################
pass_dir="./1-pass_mutect2"
mutsig_dir="./2-mutsigcv_output"
mutsig_temp_dir="$mutsig_dir/temp"
COMBINE_MAF="$mutsig_dir/combine.maf"
mkdir -p $pass_dir
mkdir -p $mutsig_dir
mkdir -p $mutsig_temp_dir
HAS_HEADER="yes"
for sample in `cat ${IDS_LIST}`;
do 
  INPUT_VCF=$VCF_PATH/1-mutect2_output/${sample}_mutect2_out_g.vcf
  PASS_VCF=$pass_dir/${sample}_pass_mutect2.vcf
  TEMP_VCF=$mutsig_temp_dir/temp_${sample}_somatic.vcf
  OUTPUT_MAF=$mutsig_temp_dir/${sample}_somatic_vep.maf
  #####################################
  # filtering mutect2 raw output
  #####################################
  CMD="vcftools --vcf $INPUT_VCF --remove-filtered-all  --recode-INFO-all --recode --stdout"
#  CMD="java -jar $HOME/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -T SelectVariants -R $HOME/resource/DNASeq/hg19/ucsc_hg19.fa -V $INPUT_VCF -o $PASS_VCF -ef "
  echo "filtering mutect2 for $sample"
  echo $CMD
  time $CMD > $PASS_VCF
  #####################################
  # Assume that we have an output file 
  # from mutect
  #####################################
  ### rm 'chr' to make it consistent to GRCh37 format
  cat $PASS_VCF|perl -lane 's/^chr// unless /^#/;print' > $TEMP_VCF
  ### make VEP annotation and convert to MAF format
  CMD="perl $HOME/apps/cancer_genomics/format_converter/vcf2maf-master/vcf2maf.pl --input-vcf $TEMP_VCF --output-maf $OUTPUT_MAF --vep-path $VEP_PATH --vep-data $VEP_DATA --ref-fasta $VEP_DATA/homo_sapiens/84_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --tumor-id ${sample}-1 --normal-id ${sample}-2"
  echo "converting vcf to maf on $sample"
  echo "[running command]: $CMD"
  time $CMD
  #rm temp/temp_${sample}*
  if [ $HAS_HEADER == 'yes' ]
  then
    # because the vcf2maf.pl would write the first line as a comment line 
    cat $OUTPUT_MAF|head -n 2|tail -n 1 > $COMBINE_MAF
    HAS_HEADER="no"
  fi
  cat $OUTPUT_MAF|perl -lane 'print if $.>2' >> $COMBINE_MAF
done

##############################################
# calling MutSigCV
##############################################

CMD="$HOME/apps/MATLAB/r2013a/bin/matlab -nodisplay < mutsig_perform.m "
echo "calling MutSigCV on combine maf"
echo "[running command]: $CMD"
time $CMD
