##### 2016-08-22
### select non-npc germline variants list
workdir="/data/home/zuoxy/data/NPC/somatic/20160519_863closing/1-somatic_mutation/2-germline_VCF"
variant0="CombineGVCFs_1810samples_378FHNPC-32EONPC-464SPNPC-489FUNDUS-208HC-72T1D-53LHT-19XZ-56SG-39HK_recal.vcf"
cohort_name="cancer-free"
sample_list="cancer-free_exome_sample.list"

###
variant1="tmp1.vcf"
variant2="tmp2.vcf"
variant3="$cohort_name.vcf"
cd $workdir

echo "get subset raw"
java -jar /data/home/zuoxy/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -R /data/home/guoym/reference/ucsc_hg19.fa -T SelectVariants --sample_file $sample_list --variant $variant0  -o $variant1 --excludeNonVariants --removeUnusedAlternates

# left align and trim variant is for indels only, to get minimal represention of indels in a cohort
## not aviable parameter --dontTrimAlleles in GATK3.4.0
#echo "left align and trim variants"
java -jar /data/home/zuoxy/apps/GenomeAnalysisTK/3.4.0/GenomeAnalysisTK.jar -R /data/home/guoym/reference/ucsc_hg19.fa -T LeftAlignAndTrimVariants --variant $variant1  -o $variant2 --splitMultiallelics #--dontTrimAlleles

echo "QC genotypes with minDP > 5 and minGQ > 10"
vcftools --vcf $variant2 --minDP 5 --maxDP 2000 --minGQ 10 --recode --stdout > $variant3


#rm $variant{1,2}
