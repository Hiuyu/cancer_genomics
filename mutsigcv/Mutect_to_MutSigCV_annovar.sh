#!/bin/sh
module load java/1.8.0_60
##############################################
# Title: MutSigCV analysis based on multiple
# 	single vcfs from Mutect
# 	Annovar version.
# Author: Xiao-yu Zuo
# Date: 2016-06-15
# 2016-06-22:
#	A full version draft.
##############################################

WORK_PATH="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene"
IDS_LIST=$1 # please input the full path of the IDS file.
#REGION="/data/home/zuoxy/resource/DNASeq/hg19/Captured-regions_Agilent_HG19_AllExonV5_UTR.bed" #only include variants in these regions for mutsig analysis
REGION="/data/home/zuoxy/resource/DNASeq/hg19/Captured-regions_Agilent_HG19_AllExonV6_UTR.bed" #only include variants in these regions for mutsig analysis
SCRIPT_HOME="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/cancer_genomics/mutsigcv"
prepare_maf_script_pl="${SCRIPT_HOME}/prepare_maf_for_mutsig.pl"
hard_filter_script_R="${SCRIPT_HOME}/hard_filter_table.R"
##############################################
# 1. annotate VCF by using annovar
##############################################
pass_dir="$WORK_PATH/1-pass_mutect2"
mutsig_dir="$WORK_PATH/3-mutsigcv_annotate"
mutsig_temp_dir="$WORK_PATH/3-mutsigcv_annotate/temp"
annovar_dir="$WORK_PATH/2-annovar_output"
COMBINE_MAF="$mutsig_dir/combine_annovar_V5.maf"
COMBINE_KEEP_MAF="$mutsig_dir/combine_annovar_V5.keep.maf"
MUTSIG_OUTPUT_PREFIX="$mutsig_dir/mutsig_V5"

mkdir -p $pass_dir
mkdir -p $mutsig_dir
mkdir -p $mutsig_temp_dir
mkdir -p $annovar_dir
HAS_HEADER="no"

for sample in `cat ${IDS_LIST}`;
do
  INPUT_VCF="$WORK_PATH/1-mutect2_output/${sample}_mutect2_out_g.vcf"
  PASS_VCF="$pass_dir/${sample}_pass_mutect2.vcf"
  ANV_OUT="$annovar_dir/$sample"
  TARGET_TABLE="$mutsig_temp_dir/${sample}.target.txt"
  #####################################
  # filtering mutect2 raw output
  # restricted to targeted region
  #####################################
  CMD="vcftools --vcf $INPUT_VCF  --remove-filtered-all  --recode-INFO-all --recode --stdout"
## seems to slow and give up 
## CMD="java -jar $HOME/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -T SelectVariants -R $HOME/resource/DNASeq/hg19/ucsc_hg19.fa -V $INPUT_VCF -o $PASS_VCF -ef "
  #echo "filtering mutect2 for $sample"
  #echo $CMD
  #time $CMD > $PASS_VCF
  #####################################
  # Assume that we have an output file
  # from mutect
  #####################################
  ### make annovar annotation and save to table format
  CMD="perl /data/home/zuoxy/apps/annovar/2015Dec14/table_annovar.pl $PASS_VCF /data/home/zuoxy/resource/annovar_hg19/humandb/ -buildver hg19 -out $ANV_OUT -remove -protocol refGene,avsnp144,cosmic74,dbnsfp30a,clinvar_20150629,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_sas,exac03,esp6500siv2_all,dbnsfp31a_interpro,dbscsnv11,kaviar_20150923,hrcr1 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --onetranscript -nastring . -vcfinput --dot2underline -thread 4"
  #echo "annotating $sample"
  #echo "[running command]: $CMD"
  #time $CMD
  #########################################
  # vcf back to table, restricted to target regions !  
  #########################################
  CMD="java -jar $HOME/apps/GenomeAnalysisTK/3.5.0/GenomeAnalysisTK.jar -T VariantsToTable -R $HOME/resource/DNASeq/hg19/ucsc_hg19.fa -V ${ANV_OUT}.hg19_multianno.vcf --allowMissingData --showFiltered -o ${TARGET_TABLE} -F CHROM -F POS -F REF -F ALT -F Gene_refGene -F Func_refGene -F avsnp144 -F ExonicFunc_refGene -F AAChange_refGene -F Kaviar_AC -F Kaviar_AN -F Kaviar_AF -F HRC_AC -F HRC_AN -F HRC_AF -F 1000g2015aug_all -F ExAC_ALL -F SIFT_pred -F SIFT_score -F PROVEAN_pred -F PROVEAN_score -GF GT -GF AD -GF AF -GF QSS  -AMD -L ${REGION}"
  #echo "vcf back to the table for $sample"
  #echo "[running command]: $CMD"
  #time $CMD
  #########################################
  # recode variant type and prepare for 
  # mutsigcv input
  #########################################
  ### because we need the first line as header
  if [ $HAS_HEADER == 'no' ]
  then
    perl $prepare_maf_script_pl $sample $TARGET_TABLE 1000g2015aug_all ExAC_ALL HRC_AF Kaviar_AF TUMOR.GT TUMOR.AD TUMOR.AF NORMAL.GT NORMAL.AD NORMAL.AF  > $COMBINE_MAF
    HAS_HEADER="yes"
  else
    perl $prepare_maf_script_pl  $sample $TARGET_TABLE 1000g2015aug_all ExAC_ALL HRC_AF Kaviar_AF TUMOR.GT TUMOR.AD TUMOR.AF NORMAL.GT NORMAL.AD NORMAL.AF|perl -lane 'print if $. >1' >> $COMBINE_MAF
  fi
done

##############################################
# do some filter, 
# e.g. population frequency,
#      sequencing coverage(depth),
#      mutant allele frequency in tumor and normal
##############################################
module load R
Rscript $hard_filter_script_R $COMBINE_MAF $COMBINE_KEEP_MAF "(ExAC_ALL<0.01 & X1000g2015aug_all<0.01 & Kaviar_AF<0.01 & HRC_AF<0.01) & TUMOR.DP > 14 & TUMOR.AD>5 & TUMOR.AF>0.05 &NORMAL.DP>8 & NORMAL.AD<5 & NORMAL.AF<0.03"


##############################################
# calling MutSigCV
##############################################
CMD="$HOME/apps/MATLAB/r2013a/bin/matlab -nodesktop -nosplash -nodisplay  -nojvm -r \"mutsig_perform('$COMBINE_KEEP_MAF','$MUTSIG_OUTPUT_PREFIX')\""
echo "calling MutSigCV on combine maf"
echo "[running command]: $CMD"
echo "cd $SCRIPT_HOME" > .mutsig.run.sh
echo $CMD >> .mutsig.run.sh
#time bash .mutsig.run.sh
#rm .mutsig.run.sh 


