#!/bin/bash
##################################
# Title:combine individual mutect/mutect2 output,
#       make annotation, and convert back to table
# Author: Xiao-yu Zuo
# Date: 2016-09-23
##################################
# grep ^# 3-mutect2_PASS/WES00003_pass_mutect2.vcf|perl -lane 'next if /^##GATK/ || /^##SAMPLE/; if(/^#CHROM/){print "##INFO=<ID=SAMPLE,Number=.,Type=String,Description=\"Sample ID\">";} print $_; ' > 3.1-mutect2_PASS/header.vcf
#########################
# Step 1: combine mutect/mutect2 output to a big vcf
#         remove chrM and other small contigs
# Assume the header file is copied here.
#########################
workdir="$PWD" # working path, all required files should be manually copied here
source_single_path="$PWD/../3-mutect2_PASS" # 6-mutect_PASS # the path of individual pass mutect2/mutect file.
single_surfix="pass_mutect2.vcf" #"mutect_PASS.vcf" # the surfix of individual pass vcf, used to get the ${sample}_$single_surfix, e.g. WES00014_mutect_PASS.vcf
OUT_PREFIX="mutect2_pass"
PASS_VCF="$OUT_PREFIX.vcf"
sample_list="all_list.txt" # sample list to include
header_file="mutect2_header.txt" # header file, use different file for mutect and mutect2

### start
### Combine all file into a big one
cd $workdir
echo "combining vcf"
cp $header_file $PASS_VCF
test -f body_${OUT_PREFIX}.vcf && rm body_${OUT_PREFIX}.vcf
# if not append SAMPLE INFO, just add it 
for sample in `cat $sample_list`
do 
  cat $source_single_path/${sample}_$single_surfix |perl -lwane '
    next if /^#/;
    $F[7]="$F[7];SAMPLE='$sample'" unless $F[7] =~ /SAMPLE=/;
    print join "\t",@F;
  ' >> body_${OUT_PREFIX}.vcf
done
# sort and remove chrM and/or other contigs(if have?)
cat body_${OUT_PREFIX}.vcf|grep -v chrM |sort -V -k 1,1 -k 2,2n >> $PASS_VCF

##########################
# Step 2: make annotation on the big vcf
# 
##########################
ANV_PREFIX=$OUT_PREFIX
ANV_TABLE="${ANV_PREFIX}_annotate.tsv"
fasta="$HOME/resource/DNASeq/hg19/ucsc_hg19.fa"
nt=8

echo "making annotation"
perl /data/home/zuoxy/apps/annovar/2015Dec14/table_annovar.pl $PASS_VCF /data/home/zuoxy/resource/annovar_hg19/humandb/ -buildver hg19 -out $ANV_PREFIX -remove -protocol refGene,avsnp144,cosmic74,dbnsfp30a,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_sas,exac03,esp6500siv2_all,dbnsfp31a_interpro,dbscsnv11,kaviar_20150923,hrcr1,sysucc_cancer_free_841sample_20161027 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput --dot2underline -thread $nt

##########################
# Step 3: convert back to table
# 
##########################
echo "to table"
INCLUDED_FLAGS="-GF GT -GF AD -GF FA" # "-GF GT -GF AD -GF AF" for mutect2, "-GF GT -GF AD -GF FA"
java -jar $HOME/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -T VariantsToTable -R $fasta -V $ANV_PREFIX.hg19_multianno.vcf --allowMissingData --showFiltered -o $ANV_TABLE -F SAMPLE -F CHROM -F POS -F REF -F ALT -F Func_refGene -F Gene_refGene -F GeneDetail_refGene -F ExonicFunc_refGene -F AAChange_refGene -F avsnp144 -F cosmic74 -F SIFT_score -F SIFT_pred -F Polyphen2_HDIV_score -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_score -F Polyphen2_HVAR_pred -F LRT_score -F LRT_pred -F MutationTaster_score -F MutationTaster_pred -F MutationAssessor_score -F MutationAssessor_pred -F FATHMM_score -F FATHMM_pred -F PROVEAN_score -F PROVEAN_pred -F VEST3_score -F CADD_raw -F CADD_phred -F DANN_score -F fathmm-MKL_coding_score -F fathmm-MKL_coding_pred -F MetaSVM_score -F MetaSVM_pred -F MetaLR_score -F MetaLR_pred -F integrated_fitCons_score -F integrated_confidence_value -F GERP++_RS -F phyloP7way_vertebrate -F phyloP20way_mammalian -F phastCons7way_vertebrate -F phastCons20way_mammalian -F SiPhy_29way_logOdds -F 1000g2015aug_all -F 1000g2015aug_afr -F 1000g2015aug_eas -F 1000g2015aug_eur -F 1000g2015aug_amr -F 1000g2015aug_sas -F ExAC_ALL -F ExAC_AFR -F ExAC_AMR -F ExAC_EAS -F ExAC_FIN -F ExAC_NFE -F ExAC_OTH -F ExAC_SAS -F esp6500siv2_all -F Interpro_domain -F dbscSNV_ADA_SCORE -F dbscSNV_RF_SCORE -F Kaviar_AF -F Kaviar_AC -F Kaviar_AN -F HRC_AF -F HRC_AC -F HRC_AN -F HRC_non1000G_AF -F HRC_non1000G_AC -F HRC_non1000G_AN -F SYSUCC_CANCER_FREE_AC -F SYSUCC_CANCER_FREE_AF -F SYSUCC_CANCER_FREE_NC -F SYSUCC_CANCER_FREE_NF -F SYSUCC_CANCER_FREE_COVF $INCLUDED_FLAGS -AMD
