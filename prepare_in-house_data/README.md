# prepare In-house data for cancer_genomics
Hiuyu@sysucc_2016-10-26

In house data are important for evaluating the false-positive variants and mutations that are mainly resulted from the systematic bias. Possible biases are, e.g. 

* sequencing platform or library(probes)
* alignment ambiguity on low complex, highly homolous or higly polymorphic regions
* possible analysis-pipeline-specific bias

In particular, in cancer genomics where the difference of tumor/normal pair is called **somatic** mutation, false-positive mutations can largely increase the difficulty of identifying truely and biologically meaningful **driver** mutation/gene of tumorgenisis. One useful way to overcome this problem is to use a large set of **cancer-free** samples that are undergoing the exactly same process including library building, sequencing, and bioinformatic analysis. The logic behind this is if the somatic mutation are non-disease-relevant or even false-positive, they are probably not rare in a set of non-disease-of-interested samples. Like the logic of using 1000 genomes and ExAC to filtering germline variants, this in-house data can also be used to filter germline risk.
