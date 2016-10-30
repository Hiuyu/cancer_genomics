# NPC genomics project scripts repository

This project can be seperated as the following parts:

1. Detect somatic mutation and indels
2. Find driver gene and key pathway 
3. mutation signatures and their associated clinical and genetic variables
4. clinical subtype classification
5. anymore?


## 1. Detect somatic mutation and indels
1. 2016/10/30 updated: updated gemline risk filtering by using new calculation of **SYSUCC_CANCER_FREE_NF** by considering the only covered samples.
	1. variants were filtered by repeatmasker.
	2. All variants were shared by at least 2 callers, including mutect, mutect2 and varscan2
	3. variants were filterd by germline risk if their maximal population frequency across any subgroups in 1000Genome and ExAC0.3 were higher than **1%**, if they were found in > **1%** of the covered in-house cancer-free cohorts. (`NOTE: the population frequency annotation were performed by annovar, in which the ExAC data were only a subset of the entire ExAC project.`)
	4. Calculation of the mutant allele fraction(AF) is not an straighforward work, because of the differernce of different callers in estimating the mutations and counting the mutant reads or bases. But after inspected, I found only subtle difference between them, in particular between the mutect and varscan2. Therefore, I planed to use a linear model by treating mutect or mutect2 as dependent variable and others as predictor. The fitted values were then used as the estimated AF for the variants.

## 2. Find driver gene and key pathway
### Driver gene prioritization
1. We chose to use **MutsigCV** algorithm to predict the driver gene for cancer. In this analysis, we performed two schemas by using two sets of variants (both with SNV and INDEL):
	1. Included coding + UTR and intron variants, i.e. all variants within a gene.
	2. Included only coding + UTR variants, i.e. only exon variants.
	we expected that the top significant genes were consistent. MutsigCV1.2 analyses were perfomed on the online server [genepatten]("https://genepattern.broadinstitute.org/gp/pages/index.jsf")


2. We planed to use another algorithm to predict driver, the **OncodriveFM**, which has the different basic logic and consideration from the MutsigCV. As the MutSigCV was called `mutation-rate based` estimator, the OncodriveFM was called `function-impacted based` detector. The OncodriveFM was a part of the [IntoGen integrative analysis platform]("http://www.intogen.org/analysis")(Version 2.4.1). 


3. We also performed an analysis called **OncodriveMUT**, which was implemented in the [**Cancer Genomics Interpreter**]("https://www.cancergenomeinterpreter.org") on the IntoGen server. The **OncodriveMUT** algorithm was recently developed to predict whether a certain single variants was a driver mutation and to prioritize the importance of the mutations. However, no related published paper was found for detailed description of the algorithm. The schema of **OncodriveMUT** provided by the website was as follow:
	
	![schema for OncodriveMUT](https://www.cancergenomeinterpreter.org/img/oncodriveMUT_schema.png)

4. As a futher filtering of uninterested gene, we resorted to use RNAseq expression data for 21 NPC tumor samples and 10 inflammation samples. We assumed that if a gene is not expressed (or only litte expression) across tumor and inflammation, it may be less important in the tumorgenesis and therefore less interested in this project. In this step, we filtered genes with the mean **Deseq** estimated expression levels < **50** for both tumor and inflammation samples. (For individual-level expression data, see the `RNAseq_21NPCtumors_10inflam_Deseq_hanbw_cp_20161030.RData`)
	Scripts to prepare these set of genes:
	```R
	library(dplyr)
	# load RNAseq expression
	load("RNAseq_21NPCtumors_10inflam_Deseq_hanbw_cp_20161030.RData")
	# to loose the criteria, only genes with both mean and median expression lower than 50 were chosen
	notExpr.gene.tumor <- mean.RNAseq.expr.deseq.tumor %>% filter(mean < 50 & median <50) %>% select(gene_symbol)
	notExpr.gene.inflammation <- mean.RNAseq.expr.deseq.inflatmmation %>% filter(mean < 50 & median <50) %>% select(gene_symbol)
	notExpr.gene.list <- notExpr.gene.tumor$gene_symbol[notExpr.gene.tumor$gene_symbol %in% notExpr.gene.inflammation$gene_symbol]
	# save it out
	save(notExpr.gene.list, file="RNAseq_tumor_inflammation_bothNotExpr_Deseq_hanbw_cp_20161030.RData")
	```
	11661/27211 genes were included in the filtered list.

### Key pathway identification


## 3. mutation signatures and their associated clinical and genetic variables

## 4. clinical subtype classification

## 5. anymore? 
e.g. ebv analysis?


