#!/bin/sh
module load java/1.8.0_60
##############################################
# Title: MutSigCV analysis based on multiple
# 	single vcfs from Mutect
# 	Annovar version.
# Author: Xiao-yu Zuo
# Date: 2016-06-15
##############################################

VCF_PATH="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene"
IDS_LIST=$1 # please input the full path of the IDS file.
##############################################
# 1. annotate VCF by using annovar
##############################################
pass_dir="./1-pass_mutect2"
mutsig_dir="./3-mutsigcv_annotate"
annovar_dir="./2-annovar_output"
COMBINE_MAF="$mutsig_dir/combine.maf"
mkdir -p $pass_dir
#mkdir -p $mutsig_dir
mkdir -p $annovar_dir
HAS_HEADER="no"
for sample in `cat ${IDS_LIST}`;
do
  INPUT_VCF=$VCF_PATH/1-mutect2_output/${sample}_mutect2_out_g.vcf
  PASS_VCF=$pass_dir/${sample}_pass_mutect2.vcf
  ANV_OUT=$annovar_dir/$sample
  #####################################
  # filtering mutect2 raw output
  #####################################
  CMD="vcftools --vcf $INPUT_VCF --remove-filtered-all  --recode-INFO-all --recode --stdout"
#  CMD="java -jar $HOME/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -T SelectVariants -R $HOME/resource/DNASeq/hg19/ucsc_hg19.fa -V $INPUT_VCF -o $PASS_VCF -ef "
  echo "filtering mutect2 for $sample"
  echo $CMD
  #time $CMD > $PASS_VCF
  #####################################
  # Assume that we have an output file
  # from mutect
  #####################################
  ### make annovar annotation and save to table format
  CMD="perl /data/home/zuoxy/apps/annovar/2015Dec14/table_annovar.pl $PASS_VCF /data/home/zuoxy/resource/annovar_hg19/humandb/ -buildver hg19 -out $ANV_OUT -remove -protocol refGene,avsnp144,cosmic74,dbnsfp30a,clinvar_20150629,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_sas,exac03,esp6500siv2_all,dbnsfp31a_interpro,dbscsnv11,kaviar_20150923,hrcr1 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --onetranscript -nastring . -vcfinput --dot2underline -thread 4"
  echo "annotating $sample"
  echo "[running command]: $CMD"
  time $CMD
  if [ $HAS_HEADER == 'yes' ]
  then
    # because the vcf2maf.pl would write the first line as a comment line
    cat $OUTPUT_MAF|head -n 2|tail -n 1 > $COMBINE_MAF
    HAS_HEADER="no"
  fi
  #cat $OUTPUT_MAF|perl -lane 'print if $.>2' >> $COMBINE_MAF
done

##############################################
# calling MutSigCV
##############################################
exit
CMD="$HOME/apps/MATLAB/r2013a/bin/matlab -nodisplay < mutsig_perform.m "
echo "calling MutSigCV on combine maf"
echo "[running command]: $CMD"
#time $CMD
