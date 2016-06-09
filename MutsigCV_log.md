https://stackedit.io/editor#

### mutation annotation format
#### headers
|column|name|
|---:|:---|
| 1| Hugo_Symbol|
| 2| Entrez_Gene_Id|
| 3| Center|
| 4| NCBI_Build|
| 5| Chromosome|
| 6| Start_Position|
| 7| End_Position|
| 8| Strand|
| 9| Variant_Classification|
| 10| Variant_Type|
| 11| Reference_Allele|
| 12| Tumor_Seq_Allele1|
| 13| Tumor_Seq_Allele2|
| 14| dbSNP_RS|
| 15| dbSNP_Val_Status|
| 16| Tumor_Sample_Barcode|
| 17| Matched_Norm_Sample_Barcode|
| 18| Match_Norm_Seq_Allele1|
| 19| Match_Norm_Seq_Allele2|
| 20| Tumor_Validation_Allele1|
| 21| Tumor_Validation_Allele2|
| 22| Match_Norm_Validation_Allele1|
| 23| Match_Norm_Validation_Allele2|
| 24| Verification_Status|
| 25| Validation_Status|
| 26| Mutation_Status|
| 27| Sequencing_Phase|
| 28| Sequence_Source|
| 29| Validation_Method|
| 30| Score|
| 31| BAM_File|
| 32| Sequencer|
| 33| Tumor_Sample_UUID|
| 34| Matched_Norm_Sample_UUID|
| 35| HGVSc|
| 36| HGVSp|
| 37| HGVSp_Short|
| 38| Transcript_ID|
| 39| Exon_Number|
| 40| t_depth|
| 41| t_ref_count|
| 42| t_alt_count|
| 43| n_depth|
| 44| n_ref_count|
| 45| n_alt_count|
| 46| all_effects|
| 47| Allele|
| 48| Gene|
| 49| Feature|
| 50| Feature_type|
| 51| Consequence|
| 52| cDNA_position|
| 53| CDS_position|
| 54| Protein_position|
| 55| Amino_acids|
| 56| Codons|
| 57| Existing_variation|
| 58| ALLELE_NUM|
| 59| DISTANCE|
| 60| STRAND|
| 61| SYMBOL|
| 62| SYMBOL_SOURCE|
| 63| HGNC_ID|
| 64| BIOTYPE|
| 65| CANONICAL|
| 66| CCDS|
| 67| ENSP|
| 68| SWISSPROT|
| 69| TREMBL|
| 70| UNIPARC|
| 71| RefSeq|
| 72| SIFT|
| 73| PolyPhen|
| 74| EXON|
| 75| INTRON|
| 76| DOMAINS|
| 77| GMAF|
| 78| AFR_MAF|
| 79| AMR_MAF|
| 80| ASN_MAF|
| 81| EAS_MAF|
| 82| EUR_MAF|
| 83| SAS_MAF|
| 84| AA_MAF|
| 85| EA_MAF|
| 86| CLIN_SIG|
| 87| SOMATIC|
| 88| PUBMED|
| 89| MOTIF_NAME|
| 90| MOTIF_POS|
| 91| HIGH_INF_POS|
| 92| MOTIF_SCORE_CHANGE|
| 93| IMPACT|
| 94| PICK|
| 95| VARIANT_CLASS|
| 96| TSL|
| 97| HGVS_OFFSET|
| 98| PHENO|
| 99| MINIMISED|
| 100| ExAC_AF|
| 101| ExAC_AF_AFR|
| 102| ExAC_AF_AMR|
| 103| ExAC_AF_EAS|
| 104| ExAC_AF_FIN|
| 105| ExAC_AF_NFE|
| 106| ExAC_AF_OTH|
| 107| ExAC_AF_SAS|
| 108| GENE_PHENO|
| 109| FILTER|

### mutsig output mutation file
```bash
[zuoxy@c035 2-mutsigcv_output]$ head -n 1  mutsig_output.mutations.txt |tr '\t' '\n'|cat -n |perl -lane 's/\s+/|/g;print$_,"|"'
```
|column|name|
|---:|:---|
|**1**|**Hugo_Symbol**|
|2|Entrez_Gene_Id|
|3|Center|
|4|NCBI_Build|
|5|Chromosome|
|6|Start_Position|
|7|End_Position|
|8|Strand|
|**9**|**Variant_Classification**|
|**10**|**Variant_Type**|
|11|Reference_Allele|
|12|Tumor_Seq_Allele1|
|13|Tumor_Seq_Allele2|
|14|dbSNP_RS|
|15|dbSNP_Val_Status|
|16|Tumor_Sample_Barcode|
|17|Matched_Norm_Sample_Barcode|
|18|Match_Norm_Seq_Allele1|
|19|Match_Norm_Seq_Allele2|
|20|Tumor_Validation_Allele1|
|21|Tumor_Validation_Allele2|
|22|Match_Norm_Validation_Allele1|
|23|Match_Norm_Validation_Allele2|
|24|Verification_Status|
|25|Validation_Status|
|26|Mutation_Status|
|27|Sequencing_Phase|
|28|Sequence_Source|
|29|Validation_Method|
|30|Score|
|31|BAM_File|
|32|Sequencer|
|33|Tumor_Sample_UUID|
|34|Matched_Norm_Sample_UUID|
|35|HGVSc|
|36|HGVSp|
|37|HGVSp_Short|
|38|Transcript_ID|
|39|Exon_Number|
|40|t_depth|
|41|t_ref_count|
|42|t_alt_count|
|43|n_depth|
|44|n_ref_count|
|45|n_alt_count|
|46|all_effects|
|47|Allele|
|48|Gene|
|49|Feature|
|50|Feature_type|
|51|Consequence|
|52|cDNA_position|
|53|CDS_position|
|54|Protein_position|
|55|Amino_acids|
|56|Codons|
|57|Existing_variation|
|58|ALLELE_NUM|
|59|DISTANCE|
|60|STRAND|
|61|SYMBOL|
|62|SYMBOL_SOURCE|
|63|HGNC_ID|
|64|BIOTYPE|
|65|CANONICAL|
|66|CCDS|
|67|ENSP|
|68|SWISSPROT|
|69|TREMBL|
|70|UNIPARC|
|71|RefSeq|
|72|SIFT|
|73|PolyPhen|
|74|EXON|
|75|INTRON|
|76|DOMAINS|
|77|GMAF|
|78|AFR_MAF|
|79|AMR_MAF|
|80|ASN_MAF|
|81|EAS_MAF|
|82|EUR_MAF|
|83|SAS_MAF|
|84|AA_MAF|
|85|EA_MAF|
|86|CLIN_SIG|
|87|SOMATIC|
|88|PUBMED|
|89|MOTIF_NAME|
|90|MOTIF_POS|
|91|HIGH_INF_POS|
|92|MOTIF_SCORE_CHANGE|
|93|IMPACT|
|94|PICK|
|95|VARIANT_CLASS|
|96|TSL|
|97|HGVS_OFFSET|
|98|PHENO|
|99|MINIMISED|
|100|ExAC_AF|
|101|ExAC_AF_AFR|
|102|ExAC_AF_AMR|
|103|ExAC_AF_EAS|
|104|ExAC_AF_FIN|
|105|ExAC_AF_NFE|
|106|ExAC_AF_OTH|
|107|ExAC_AF_SAS|
|108|GENE_PHENO|
|109|FILTER|
|**110**|**gene**|
|111|patient|
|**112**|**effect**|
|113|chr|
|114|start|
|115|ref_allele|
|**116**|**categ**|

### variant\_classification vs mutation\_type_dictionary\_file
#### *variant_classification* column in **maf** format
>|count|variant_classification
|---|---
| 5155| 3'Flank
| 35849| 3'UTR
| 2470| 5'Flank
| 4794| 5'UTR
| 1296| Frame_Shift_Del
| 441| Frame_Shift_Ins
| 295| IGR
| 194| In_Frame_Del
| 215| In_Frame_Ins
| 36641| Intron
| 27841| Missense_Mutation
| 1074| Nonsense_Mutation
| 20| Nonstop_Mutation
| 2344| RNA
| 20494| Silent
| 331| Splice_Site
| 3455| Targeted_Region
| 41| Translation_Start_Site

Following is an explanation of the variant_classification from [VEP](!"http://asia.ensembl.org/info/genome/variation/predicted_data.html#consequences") and [this post](!"https://www.biostars.org/p/69222/").
> **how to understand the variant_classification* in MAF**
> See below a diagram showing the location of each display term relative to the transcript structure
![a diagram showing the location of each display term relative to the transcript structure](http://asia.ensembl.org/info/genome/variation/consequences.jpg)
> 

#### *mutation\_type_dictionary\_file* in **MutSigCV** dependencies folder

Variant_Classification|effect
---|---
Silent|silent
Synonymous|silent
Missense|nonsilent
Missense_Mutation|nonsilent
Nonsense|null
Nonsense_Mutation|null
Nonstop_Mutation|null
Read-through|null
Frame_Shift_Del|null
Frame_Shift_Ins|null
In_Frame_Del|null
In_Frame_Ins|null
Splice|null
Splice_Region|null
Splice_Site|null
Splice_Site_Del|null
Splice_Site_DNP|null
Splice_Site_Ins|null
Splice_Site_ONP|null
Splice_Site_SNP|null
Start_Codon_Del|null
Start_Codon_DNP|null
Start_Codon_Ins|null
Start_Codon_ONP|null
Stop_Codon_Del|null
Stop_Codon_DNP|null
Stop_Codon_Ins|null
Translation_Start_Site|null
De_novo_Start|null
De_novo_Start_InFrame|null
De_novo_Start_OutOfFrame|null
IGR|noncoding
Intron|noncoding
3'Flank|noncoding
3'Promoter|noncoding
3'UTR|noncoding
3'-UTR|noncoding
5'Flank|noncoding
5'-Flank|noncoding
5'Promoter|noncoding
5'UTR|noncoding
5'-UTR|noncoding
downstream|noncoding
miRNA|noncoding
NCSD|noncoding
Non-coding_Transcript|noncoding
Promoter|noncoding
RNA|noncoding
upstream|noncoding
upstream;downstream|noncoding

### Fnding some differrence
#### Get information from The `mutation.txt` generated by MutsigCV
```bash
[zuoxy@c035 2-mutsigcv_output]$ pwd
/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-mutsigcv_output
[zuoxy@c035 2-mutsigcv_output]$ cut -f 1,9,10,11,12,13,40-45,51,56,72-75,110,112,115,116 mutsig_output.mutations.txt > minor.mut.txt
[zuoxy@c035 2-mutsigcv_output]$ head minor.mut.txt |perl -lane '@F=split /\t/;print "---|"x@F if $.==2;print join "|",@F;'
```
Hugo_Symbol|Variant_Classification|Variant_Type|Reference_Allele|Tumor_Seq_Allele1|Tumor_Seq_Allele2|t_depth|t_ref_count|t_alt_count|n_depth|n_ref_count|n_alt_count|Consequence|Codons|SIFT|PolyPhen|EXON|INTRON|gene|effect|ref_allele|categ
---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
MST1L|RNA|SNP|G|G|A|60|56|4|42|41|1|non_coding_transcript_exon_variant||||10/20||MST1L|noncoding|G|*CpG->T
CYP4A22|Intron|SNP|C|C|T|209|204|5|155|155|0|intron_variant|||||11/11|CYP4A22|noncoding|C|*Np(A/C/T)->transit
LYSMD1|5'UTR|SNP|G|G|A|214|208|6|110|110|0|5_prime_UTR_variant||||1/3||LYSMD1|noncoding|G|*CpG->T
SMG7|Missense_Mutation|SNP|G|G|C|152|130|22|94|94|0|missense_variant|Gtt/Ctt|tolerated_low_confidence(0.23)|benign(0.328)|17/23||SMG7|nonsilent|G|transver
SRP9|3'UTR|SNP|G|G|A|231|224|7|158|157|1|3_prime_UTR_variant||||3/3||SRP9|noncoding|G|*Np(A/C/T)->transit
OBSCN|Silent|SNP|G|G|A|124|116|8|52|52|0|synonymous_variant|caG/caA|||51/116||OBSCN|silent|G|*Np(A/C/T)->transit
HS6ST1|Missense_Mutation|SNP|G|G|T|24|17|7|12|12|0|missense_variant|Cgc/Agc|deleterious(0)|probably_damaging(0.998)|2/2||HS6ST1|nonsilent|G|transver
CCDC74B|Missense_Mutation|SNP|G|G|A|38|33|5|10|10|0|missense_variant|cCg/cTg|tolerated(0.23)|benign(0.028)|6/8||CCDC74B|nonsilent|G|*CpG->T
SMPD4|3'UTR|SNP|G|G|T|99|95|4|86|85|1|3_prime_UTR_variant||||20/20||SMPD4|noncoding|G|transver

#### Get information from The `mutation.txt` generated by MutsigCV
```bash
[zuoxy@c035 2-mutsigcv_output]$ pwd
/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-mutsigcv_output
[zuoxy@c035 2-mutsigcv_output]$ cut -f 1,9,10,11,12,13,40-45,51,56,72-75,110,112,115,116 mutsig_output.mutations.txt > minor.mut.txt
[zuoxy@c035 2-mutsigcv_output]$ head minor.mut.txt |perl -lane '@F=split /\t/;print "---|"x@F if $.==2;print join "|",@F;'
```
