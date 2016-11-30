library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
#library(VennDiagram)
library(GenomicRanges)
# in windows
#setwd("C:/Users/LibJu/workspace/myData/NPC/cancer_genomics/caller_comparison")
# setwd("C:/Users/LibJu/workspace/myData/P1-NPC/cancer_genomics/1-somatic_mutation/temp/")
#source("C:/Users/LibJu/workspace/Github-project_scripts-repository/cancer_genomics/NPC_genomics/0-helper_function.R")
# in cluster
source("/data/home/zuoxy/data/NPC/somatic/20160519_863closing/11-oncomap_analysis_workspace/scripts_npc_project_also_in_github/cancer_genomics/NPC_genomics/0-helper_function.R")

####
#start
####
message("start")
mutect2<-as.data.frame(data.table::fread("mutect2_pass_annotate.tsv"))
mutect<-as.data.frame(data.table::fread("mutect_pass_annotate.tsv"))
varscan2<-as.data.frame(data.table::fread("varscan2_pass_annotate.tsv"))
## get the population frequency colnames
pop.freq.col=colnames(mutect2)[sapply(colnames(mutect2),function(x){
  grepl("ExAC",x)|grepl("1000g",x)|grepl("CANCER_FREE_[AN]F",x)
})]

# if no annotation for the Func_refGene and Gene_refGene, remove the variant
mutect2 = mutect2 %>% filter(Func_refGene != "." & Gene_refGene != ".")
mutect = mutect %>% filter(Func_refGene != "." & Gene_refGene != ".")
varscan2 = varscan2 %>% filter(Func_refGene != "." & Gene_refGene != ".")

print("removing ambigious annotation")
amb=with(mutect2, wrapper_remove_ambigious_find_nearest(Gene_refGene,Func_refGene,ExonicFunc_refGene,GeneDetail_refGene,sep1="\\\\x3b",sep2="\\\\x3d"))
mutect2=within(mutect2,{
  source=decide_source(SAMPLE)
  type=decide_variant_type(REF,ALT)
  Gene_refGene=amb$gene
  Func_refGene=amb$func_gene
  ExonicFunc_refGene=amb$func_exonic
  GeneDetail_refGene=amb$dist
})
amb=with(mutect, wrapper_remove_ambigious_find_nearest(Gene_refGene,Func_refGene,ExonicFunc_refGene,GeneDetail_refGene,sep1="\\\\x3b",sep2="\\\\x3d"))
mutect=within(mutect,{
  source=decide_source(SAMPLE)
  type=decide_variant_type(REF,ALT)
  Gene_refGene=amb$gene
  Func_refGene=amb$func_gene
  ExonicFunc_refGene=amb$func_exonic
  GeneDetail_refGene=amb$dist
})
varscan2$ExonicFunc_refGene=gsub("\\s+","_",varscan2$ExonicFunc_refGene)
amb=with(varscan2, wrapper_remove_ambigious_find_nearest(Gene_refGene,Func_refGene,ExonicFunc_refGene,GeneDetail_refGene,sep1=";",sep2="="))
varscan2=within(varscan2,{
  source=decide_source(SAMPLE)
  type=decide_variant_type(REF,ALT)
  Gene_refGene=amb$gene
  Func_refGene=amb$func_gene
  ExonicFunc_refGene=amb$func_exonic
  GeneDetail_refGene=amb$dist
})
## if for HK and SG, the source field was set as follows:
# > mutect$source=gsub("^(\\w+)_.+","\\1",mutect$SAMPLE)
# > mutect2$source=gsub("^(\\w+)_.+","\\1",mutect2$SAMPLE)
# > varscan2$source=gsub("^(\\w+)_.+","\\1",varscan2$SAMPLE)
save.image("raw_all_input_3caller_3callers.RData")
####### part 1 WES SNV comparison #####
# remove some discordant columns to make them "rbind"-able
D.mt2=mutect2%>%select(everything(),-TUMOR.FA,-NORMAL.FA)
D.mt=mutect%>%select(everything(),-TUMOR.FA,-NORMAL.FA)
D.vs=varscan2%>%select(everything(),POS=Start,-End)

# total number of SNV count
print(paste("ALL variant count for mutect2, mutect, varscan is ",nrow(D.mt2), nrow(D.mt), nrow(D.vs)))
caller=c(rep("mutect2",nrow(D.mt2)),
         rep("mutect",nrow(D.mt)),
         rep("varscan2",nrow(D.vs))
)
D.all=cbind(rbind(D.mt2,D.mt,D.vs),caller)

## recode population frequency to be numeric
for(x in pop.freq.col){
  D.all[D.all[,x]==".",x]="0"
  D.all[,x]=as.numeric(D.all[,x])
  D.all[is.na(D.all[,x]),x]=0
}

# very time consuming ! ! ! but no need in cluster
# find the max pop.freq in 1000G and ExAC
# tmp.freq=D.all[,pop.freq.col]
# D.all$max.pop.freq=0
# print("finding max pop freq")
# i=1
# while(i<=nrow(D.all)){
#   D.all$max.pop.freq[i]=max(tmp.freq[i,])
#   if(i%%50000==0){message(paste0("parsed ",i/1000,"K sites"))}
#   i=i+1
# }
print("finding max pop freq")
tmp.freq=apply(D.all[,pop.freq.col],1,function(x) max(x))
D.all$max.pop.freq = tmp.freq

D.all$pop.freq=.assign_pop_freq_group(D.all$max.pop.freq)
D.all_snv_indel=D.all
# set index
D.all_snv_indel=D.all_snv_indel%>%mutate(index=paste0(SAMPLE,"_",CHROM,":",POS,"_",REF,">",ALT))

# add region column flag
D.all_snv_indel$region="non-exon"
D.all_snv_indel$region[D.all_snv_indel$Func_refGene %in% c("exonic","UTR3","UTR5","splicing")] = "exon+UTR+splicing"
D.all_snv_indel=within(D.all_snv_indel,{ExonicFunc_refGene[ExonicFunc_refGene == "."]="non-exonic"})
# temp script for plotting
#tmp.df=with(D.all_snv_indel%>%filter(ExonicFunc_refGene!="."),prop.table(table(ExonicFunc_refGene,caller,region,type,source),2:5))
#tmp.df=melt(tmp.df)
#ggplot(tmp.df, aes(x=caller,y=value,fill=ExonicFunc_refGene))+geom_bar(stat="identity",position="stack")+ylab("frequency") +
#  ggtitle("only Exonic variants") +facet_grid(type+region~source)

## saving objects
output_file="ALL_SNV_INDEL_3callers_updated_in-house-data_20161028.RData"
save(D.all_snv_indel,file=output_file)
print(paste("all done, please load()",output_file,"in R"))

# ## generate concordance table for callers
# tmp.rmcol=colnames(D.all_snv_indel)[grepl("AFR|EAS|EUR|AMR|FIN|OTH|SAS|NFE$",colnames(D.all_snv_indel),ignore.case = TRUE)]
# D.all.1=D.all_snv_indel%>%select(-one_of(tmp.rmcol),-starts_with("TUMOR"),-starts_with("NORMAL"))
# # reconstruct data
# D.all.2=subset(D.all.1,!duplicated(D.all.1$index))
# D.all.2$mutect=FALSE; D.all.2$mutect[D.all.2$index %in% D.all.1$index[D.all.1$caller=="mutect"]]=TRUE
# D.all.2$mutect2=FALSE; D.all.2$mutect2[D.all.2$index %in% D.all.1$index[D.all.1$caller=="mutect2"]]=TRUE
# D.all.2$varscan2=FALSE; D.all.2$varscan2[D.all.2$index %in% D.all.1$index[D.all.1$caller=="varscan2"]]=TRUE
# D.all.2$share_times= D.all.2$mutect + D.all.2$mutect2 + D.all.2$varscan2

# ## filter germline
# D.all.1.filter=D.all.1%>%filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3 )
# D.all.2.filter=D.all.2%>%filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3 )
#D.all.2.exonic %>% filter(Func_refGene %in% c("UTR3","UTR5","splicing","exonic"))

# ## concordance analysis: venn plot
# # all variants, before filtering
# tmp.d=D.all.2
# venn.diagram(
#   list(
#     mutect=which(with(tmp.d,mutect & type=="SNV")),
#     mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
#     varscan=which(with(tmp.d,varscan2 & type=="SNV"))
#   ),
#   main="All SNVs, before germline filtering",
#   "venn_2caller_snv_all.tiff"
# )
# venn.diagram(
#   list(
#     mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
#     varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
#   ),
#   main="All INDELs, before germline filtering",
#   "venn_2caller_indel_all.tiff"
# )
# # all variants , after filtering
# tmp.d=D.all.2
# venn.diagram(
#   list(
#     mutect=which(with(tmp.d,mutect & type=="SNV")),
#     mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
#     varscan=which(with(tmp.d,varscan2 & type=="SNV"))
#   ),
#   main="All SNVs, after germline filtering",
#   "venn_2caller_snv_filter.tiff"
# )
# venn.diagram(
#   list(
#     mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
#     varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
#   ),
#   main="All INDELs, after germline filtering",
#   "venn_2caller_indel_filter.tiff"
# )
# # exons only, before filtering
# tmp.d=D.all.2%>%filter(region=="exon+UTR+splicing")
# venn.diagram(
#   list(
#     mutect=which(with(tmp.d,mutect & type=="SNV")),
#     mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
#     varscan=which(with(tmp.d,varscan2 & type=="SNV"))
#   ),
#   main="All SNVs, before germline filtering",
#   "venn_2caller_snv_all_exon.tiff"
# )
# venn.diagram(
#   list(
#     mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
#     varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
#   ),
#   main="All INDELs, before germline filtering",
#   "venn_2caller_indel_all_exon.tiff"
# )
# # exons only , after filtering
# tmp.d=D.all.2.filter%>%filter(region=="exon+UTR+splicing")
# venn.diagram(
#   list(
#     mutect=which(with(tmp.d,mutect & type=="SNV")),
#     mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
#     varscan=which(with(tmp.d,varscan2 & type=="SNV"))
#   ),
#   main="All SNVs, after germline filtering",
#   "venn_2caller_snv_filter_exon.tiff"
# )
# venn.diagram(
#   list(
#     mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
#     varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
#   ),
#   main="All INDELs, after germline filtering",
#   "venn_2caller_indel_filter_exon.tiff"
# )
# save.image("workspace_call_comparison_20160926.RData")


## repeat region (repeatmasker) overlap
#load("/data/home/zuoxy/resource/genome/ucsc/hg19/annotation/Repeats/rmsk_GRange.RData")
load("rmsk_GRange.RData")
a=with(D.all_snv_indel,convert_0base_to_1base(POS,REF,ALT))
D.all_snv_indel=within(D.all_snv_indel,{
  start=a$START
  end=a$END
  REF1=a$NEW.REF
  ALT1=a$NEW.ALT
  strand=rep("+",nrow(D.all_snv_indel))
})
D.all_snv_indel_grange = makeGRangesFromDataFrame(D.all_snv_indel,seqnames.field="CHROM",keep.extra.columns=T,starts.in.df.are.0based = FALSE)
tmp.rpid=subsetByOverlaps(D.all_snv_indel_grange,rmsk_grange,type="any",ignore.strand=T)
D.norep_snv_indel_grange=D.all_snv_indel_grange[!D.all_snv_indel_grange$index %in% tmp.rpid$index,]

## gene-level concordance
# raw data
# tmp.d=D.all.2%>%filter(source=="WES" & Func_refGene!="intergenic")
# tmp.d1=D.all.2 %>% filter(index %in% unique(D.norep_snv_indel_grange$index)) %>%
#   filter(source=="WES" & Func_refGene!="intergenic") %>%
#   filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3)
# venn.diagram(
#   list(
#     mutect=unique(tmp.d%>%filter(mutect)%>%select(Gene_refGene))[,1],
#     mutect2=unique(tmp.d%>%filter(mutect2)%>%select(Gene_refGene))[,1],
#     varscan2=unique(tmp.d%>%filter(varscan2)%>%select(Gene_refGene))[,1]
#   ),
#   main="All genes, after region & germline filtering",
#   "venn_3caller_gene_filter.tiff"
# )

## liuwsh heatmap gene comparison
# liu.data=liu.data%>%mutate(index=paste0(SAMPLE,"_",Chr,":",Start,"-",End,"_",Ref,">",Alt))
# D.all_snv_indel=as.data.frame(D.all_snv_indel_grange)
# D.all_snv_indel=D.all_snv_indel%>%mutate(index=paste0(SAMPLE,"_",seqnames,":",start,"-",end,"_",REF1,">",ALT1))
# tmp.d=D.all.2%>%filter(index%in%liu.data$index)
# venn.diagram(
#   list(
#     mutect=unique(tmp.d%>%filter(mutect)%>%select(index))[,1],
#     mutect2=unique(tmp.d%>%filter(mutect2)%>%select(index))[,1],
#     varscan2=unique(tmp.d%>%filter(varscan2)%>%select(index))[,1]
#   ),
#   main="liu, variant",
#   "liu.venn_3caller_variant.tiff"
# )
# tmp.d_grange=makeGRangesFromDataFrame(tmp.d,keep.extra.columns = T)
# tmp.rp=subsetByOverlaps(tmp.d_grange,rmsk_grange,type="any",ignore.strand=T)


## top 50 high frequency gene.
# z1.5=D.all.2%>%filter(share_times>0 & mutect2) %>%
#   filter(source=="WES") %>% 
#   filter(region!="non-exon") %>%
#   filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>%
#   group_by(Gene_refGene) %>% summarise(n=n()) %>% arrange(-n)

## only norep variants.
D.norep_snv_indel=as.data.frame(D.norep_snv_indel_grange) %>% select(CHROM=seqnames,everything())
tmp.rmcol=colnames(D.norep_snv_indel)[grepl("AFR|EAS|EUR|AMR|FIN|OTH|SAS|NFE$",colnames(D.norep_snv_indel),ignore.case = TRUE)]
D.norep.1=D.norep_snv_indel%>%select(-one_of(tmp.rmcol),-starts_with("TUMOR"),-starts_with("NORMAL"))
D.norep.2=subset(D.norep.1,!duplicated(D.norep.1$index))
D.norep.2$mutect=FALSE; D.norep.2$mutect[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="mutect"]]=TRUE
D.norep.2$mutect2=FALSE; D.norep.2$mutect2[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="mutect2"]]=TRUE
D.norep.2$varscan2=FALSE; D.norep.2$varscan2[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="varscan2"]]=TRUE
D.norep.2$share_times= D.norep.2$mutect + D.norep.2$mutect2 + D.norep.2$varscan2
#save(D.norep.1,D.norep.2,file="norep_3callers_20160930.RData")

## exclude excluded or not-interested gene (PMID:22294350)
#backlist_gene_PMID22294350=read.table("PMID22294350_Table_S7_gene_exclusion_list_final.txt",stringsAsFactors = F,header=F,sep="\t")
#backlist_gene_PMID22294350=toupper(backlist_gene_PMID22294350[,1])
#save(backlist_gene_PMID22294350,file="backlist_gene_PMID22294350_TableS7.RData")
load("backlist_gene_PMID22294350_TableS7.RData")
# excluded backlist genes and MUC family,
# D.norep.noback.exon.2=D.norep.2%>%filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
#   filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
#   filter(! grepl("\\bOR\\d+",Gene_refGene)) %>%
#   filter(region=="exon+UTR+splicing")

# ## make gene summary for several data
# # z1: WES 3 callers; z2: WES 2 callers
# z1.summary = summarize_by_gene_annovar(D.norep.noback.exon.2 %>% filter((share_times>2 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
#                                          filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>% 
#                                          filter(source=="WES") %>% filter(Func_refGene%in%c("exonic","splicing"))
#                                        )
# z2.summary = summarize_by_gene_annovar(D.norep.noback.exon.2 %>% filter((share_times>1 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
#                                          filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>%
#                                          filter(source=="WES") %>% filter(Func_refGene%in%c("exonic","splicing"))
#                                        )
# write.csv(z1.summary,"WES_3caller.csv")
# write.csv(z2.summary,"WES_2caller.csv")
# 
# z1.1=z1.summary[,1:2];z1.1$id=1:nrow(z1.1)
# z2.1=z2.summary[,1:2];z2.1$id=1:nrow(z2.1)
# order=rep(c("3caller","2caller"),c(nrow(z1.summary),nrow(z2.summary)))
# z12=cbind(rbind(z1.1,z2.1),order)
# z12.1=dcast(z12,Gene_refGene~order,value.var = "n_sample")
# colnames(z12.1)=c("gene","n.2caller","n.3caller")
# z12.1=z12.1%>%arrange(gene)
# z12.2=dcast(z12,Gene_refGene~order,value.var = "id")
# colnames(z12.2)=c("gene","order.2caller","order.3caller")
# z12.2=z12.2%>%arrange(gene)
# z12.3=merge(z12.1,z12.2)
# z12.3=z12.3%>%arrange(order.2caller)
# write.csv(z12.3,"WES_2-3caller.csv")
# 
# tmp.z=z12.3%>%filter(order.3caller <=100 & order.2caller<=100 )
# ggplot(tmp.z,aes(order.2caller,order.3caller))+geom_point()+geom_text_repel(aes(label=gene))+
#   geom_hline(yintercept = c(50,100),lty=2)+geom_vline(xintercept = c(50,100),lty=2) +
#   geom_smooth()

# 
# library(org.Hs.eg.db)
# tmp.z12=head(z12.3,1000)
# gene <- z12.3$gene
# geneid <- select(org.Hs.eg.db, keys=gene, keytype="SYMBOL",
#                  columns="ENTREZID")



### combine all filtering steps
# single method pass QC
D.norep_snv_indel_wgs_QC=D.norep_snv_indel %>%
  filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
  filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
  filter(!grepl("\\bOR\\d+",Gene_refGene)) %>%
  filter(region=="exon+UTR+splicing") %>%
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NF < 0.01) %>% 
  filter(source=="WGS") %>%
  mutate(REF=REF1,ALT=ALT1) %>%
  select(-REF1, -ALT1,-POS) %>%
  arrange(CHROM,start)

# pass QC and shared variants
D.norep.noback.coding_UTR.wgs=D.norep.2 %>%
  filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
  filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
  filter(!grepl("\\bOR\\d+",Gene_refGene)) %>%
  filter(Func_refGene%in%c("exonic","splicing","UTR3","UTR5")) %>%
  filter((share_times>1 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NF < 0.01) %>% 
  filter(source=="WGS") %>%
  mutate(REF=REF1,ALT=ALT1) %>%
  select(-REF1, -ALT1,-POS) %>%
  arrange(CHROM,start)

## with intron variants
D.norep.noback.coding_UTR_intron.wgs=D.norep.2 %>%
  filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
  filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
  filter(!grepl("\\bOR\\d+",Gene_refGene)) %>%
  filter(Func_refGene%in%c("exonic","splicing","UTR3","UTR5","intronic")) %>%
  filter((share_times>1 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NF < 0.01) %>% 
  filter(source=="WGS") %>%
  mutate(REF=REF1,ALT=ALT1) %>%
  select(-REF1, -ALT1,-POS) %>%
  arrange(CHROM,start)

# melt form of data, with mutation AF
D.norep.noback.coding_UTR.wgs.withAD = D.norep_snv_indel %>%
  filter(index %in% D.norep.noback.coding_UTR.wgs$index)

# save all objects in QC steps
save(D.all_snv_indel, # raw variants for WGS and WES, SNV and INDEL, before any filtering
     D.norep_snv_indel, # raw variants + repeatmasker filtering
     D.norep_snv_indel_wgs_QC, # raw varaints + repeatmasker + blacklist gene + exon_UTR + germline filtering
     D.norep.noback.coding_UTR.wgs.withAD, # filtered variants + shared by 2 callers
     D.norep.noback.coding_UTR.wgs, # another form of shared variants past QC
     file="npc_norep_nobacklist_shared2_pop0.01_wgs_objects_20161028.RData"
)

##### tricks: convert annovar output to MAF. 
## require maftools packages, which is not available on cluster
# to maftools
tmp.d=D.norep.noback.coding_UTR.wgs %>%
  dplyr::select(Chr=CHROM,
                Start=start,
                End=end,
                Ref=REF,
                Alt=ALT,
                Func.refGene=Func_refGene,
                Gene.refGene=Gene_refGene,
                GeneDetail.refGene=GeneDetail_refGene,
                ExonicFunc.refGene=ExonicFunc_refGene,
                AAChange.refGene=AAChange_refGene,
                SAMPLE=SAMPLE
  )
prefix="wgs_exonic_updated_SYSUCC_NF"
out_file=paste0(prefix,".txt")
cat(colnames(tmp.d)[1:11], file=out_file,sep="\t")
cat("\n",file=out_file,append=TRUE)
cat(c(colnames(tmp.d)[1:10],"Tumor_Sample_Barcode"),file=out_file,append = TRUE,sep="\t")
cat("\n",file=out_file,append=TRUE)
tmp.d$ExonicFunc.refGene=gsub("_"," ",tmp.d$ExonicFunc.refGene)
write.table(tmp.d,out_file,sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE,append=TRUE)
###################### use maftools packages ##########################
maf_file=paste0(prefix,".maf")
var.maf=annovarToMaf(annovar = out_file, Center = 'BEI_SYSUCC', refBuild = 'hg19',
                     tsbCol = 'Tumor_Sample_Barcode', table = 'refGene', header = TRUE,MAFobj = FALSE)
write.table(var.maf,maf_file,sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
#npc=read.maf(maf = maf_file, removeSilent = FALSE, useAll = FALSE)
#tiff("bbb.tif",width=2000,height=1500)
#oncoplot(maf = npc, top = 100,removeNonMutated=FALSE,fontSize = 12)
#dev.off()
#npc.oncodrive=oncodrive(maf=npc,minMut = 5,pvalMethod = "combined")
#plotOncodrive(res = npc.oncodrive, fdrCutOff = 0.05, useFraction = FALSE)




## prepare for next step: signature and prognostics
D.coding_splicing.wgs=D.norep.noback.coding_UTR.wgs %>% 
  filter(Func_refGene %in% c("exonic", "splicing")) %>%
  filter(type=="SNV")
save(D.coding_splicing.wgs,file="WGS_SNV_coding-splicing_updated_SYSUCC_2016-10-30.RData")

## sample-gene mutation matrix
tmp.g=npc@gene.summary%>%filter(MutatedSamples>10)%>%dplyr::select(Hugo_Symbol)
tmp.g=tmp.g$Hugo_Symbol
tmp.d=D.coding_splicing.wes%>%filter(Gene_refGene%in%tmp.g)
gene.mut_mat=dcast(tmp.d,SAMPLE~Gene.refGene, value.var = "type", fun.aggregate = length)
gene.mut_mat[gene.mut_mat>1]=1
write.csv(gene.mut_mat,file = "sample_gene_matrix.csv")


load("npc_norep_nobacklist_shared2_pop0.01_wes_objects_20161028.RData")
## summary gene status, compute for all methods
## use only coding SNV
gene.summary = lapply(unique(D.norep_snv_indel_wgs_QC$caller),function(x){
  summarize_by_gene_annovar(D.norep_snv_indel_wgs_QC %>%
                              filter(caller==x)%>%
                              filter(ExonicFunc_refGene != "non-exonic" & type=="SNV") %>% 
                              filter(ExonicFunc_refGene != "unknown")
  )
})
names(gene.summary)=unique(D.norep_snv_indel_wgs_QC$caller)
gene.summary[["shared"]]=summarize_by_gene_annovar(D.norep.noback.coding_UTR.wgs %>%
                                                     filter(ExonicFunc_refGene != "non-exonic" & type=="SNV")  %>% 
                                                     filter(ExonicFunc_refGene != "unknown")
)
tmp.summary=NULL
tmp.name=names(gene.summary)
for(i in 1:length(tmp.name)){
  tmp.summary=rbind(tmp.summary,gene.summary[[tmp.name[i]]])
}
tmp.summary$method = factor(rep(names(gene.summary),sapply(gene.summary,nrow)),levels=tmp.name)

# dn/ds plot
tmp.d=tmp.summary%>%mutate(sample_frac=n_sample/474, dn_ds=(n_nonsense+n_missense)/n_silent)
tmp.max=max(tmp.d$dn_ds[!(tmp.d$dn_ds==Inf | is.nan(tmp.d$dn_ds))])
tmp.d$dn_ds[tmp.d$dn_ds==Inf]=tmp.max+5
tmp.d=tmp.d%>%filter(!is.nan(dn_ds))
intersted=c("CYLD","TP53","RYR2","TRAF3","NFKBIA","TNFAIP3","PTEN","KRAS","NRAS","ARID1A","BAP1","KMT2C","KMT2D","ATM","UBN1")
nfkb=c("TNF","RAN","TNFAIP3","CYLD","BCL10","RIPK2","PRKCA","TNFRSF1A","NFKB1","TRAF6","MALT1","ERC1","RELA","IKBKG","IKBKB","BIRC2","XPO1","CHUK","NOD2","UBE2D3","NFKBIA","ATM","BTRC")
top2=c("RHOA","DCN","EP300","TGFBR2","DGKG","PIK3C2G","PIK3CA","PLCZ1","CACNA1A","CACNA1C","DUSP6","FGF6","FGFR3","HRAS","NRAS","NTF3","TGFBR2","TP53","TRAF2","FGF23","TAOK1")
library(ggrepel)
# plot each method and shared results
plist=list()
for(s in tmp.name){
  tmp.d1 = tmp.d%>%filter(method==s)
  plist[[s]] <- ggplot(tmp.d1,aes(x=sample_frac,y=log(dn_ds,base=2)))+
    geom_point() +
    scale_x_continuous(labels = scales::percent) +
    geom_vline(xintercept = c(0.01,0.03),lty=2) +
    geom_hline(yintercept = log(c(0.25,1,4),base = 2),lty=2) +
    ggtitle(s)
}
thres=c(0.05,0.05,0.05,0.03)
plist.1=plist
plist.1=lapply(1:length(thres),function(x){
  plist[[x]] + geom_text_repel(data=tmp.d %>% filter(method==tmp.name[x]) %>%
                                 filter((sample_frac > thres[x] & (dn_ds > 4 | dn_ds <0.25)) | 
                                          Gene_refGene %in% nfkb), 
                               aes(label=Gene_refGene))
})
tmp.p=lapply(plist.1, function(x) ggplotGrob(x))
pCombine <- arrangeGrob( tmp.p[[1]],tmp.p[[2]],tmp.p[[3]],tmp.p[[4]] ,ncol=2,nrow=2)
grid.newpage()
jpeg("dn-ds_for_all_method.jpeg",width=800,height=600)
grid.draw(pCombine)
dev.off()

### the mutant allele frequency concordance between shared methods
tmp.ad = D.norep.noback.coding_UTR.wgs.withAD %>% 
  filter(ExonicFunc_refGene != "non-exonic" & type=="SNV") %>%
  select(index, caller, TUMOR.AD, NORMAL.AD)
tmp.z = parse_GATK_AD_field(tmp.ad$TUMOR.AD)
tmp.ad$TUMOR.AF=tmp.z$AF
tmp.ad$TUMOR.DP=tmp.z$DP
tmp.af = tmp.ad %>% select(everything(),-TUMOR.AD, -NORMAL.AD)
# get the AF matrix
tmp.afm = dcast(tmp.af, index~caller,value.var="TUMOR.AF")
p1=ggplot(tmp.afm,aes(mutect,mutect2))+geom_point() + ggtitle("mutect vs mutect2")
p2=ggplot(tmp.afm,aes(mutect2,varscan2))+geom_point() + ggtitle("mutect2 vs varscan2")
p3=ggplot(tmp.afm,aes(mutect,varscan2))+geom_point() + ggtitle("mutect vs varscan2")
gl <- lapply(list(p1,p2,p3), function(x) ggplotGrob(x + geom_smooth(method="lm",se=TRUE) + geom_abline(slope=1,intercept = 0) + theme_bw()))
gt <- gtable(widths=unit(rep(1,3), "null"),
             heights=unit(1, "null"))
gt <- gtable_add_grob(gt, gl, l=1:3,r=1:3,t=1,b=1)
grid.newpage()
jpeg("AF comparison_all.jpeg",width=900,height=600)
grid.draw(gt)
dev.off()

## since the AF is somewhat different for different method, using fitted value as the combined one
part1 = tmp.afm%>%filter(!(is.na(mutect)|is.na(mutect2)) & is.na(varscan2))
part2 = tmp.afm%>%filter(!(is.na(mutect)|is.na(varscan2))& is.na(mutect2))
part3 = tmp.afm%>%filter(!(is.na(mutect2)|is.na(varscan2))& is.na(mutect))
part4 = tmp.afm%>%filter(!(is.na(mutect2)|is.na(varscan2)| is.na(mutect)))
tmp.m2m=lm(mutect2~mutect,data=part4[,-1])$fitted.values
tmp.m2v2=lm(mutect2~varscan2,data=part4[,-1])$fitted.values
tmp.mv2=lm(mutect~varscan2,data=part4[,-1])$fitted.values
tmp.mm2=lm(mutect~mutect2,data=part4[,-1])$fitted.values
tmp.v2m2=lm(varscan2~mutect2,data=part4[,-1])$fitted.values
tmp.v2m=lm(varscan2~mutect,data=part4[,-1])$fitted.values
resuid=cbind(tmp.m2m,tmp.mm2,tmp.mv2,tmp.v2m,tmp.m2v2,tmp.v2m2)
colnames(resuid)=c("m2_m","m_m2","m_v2","v2_m","m2_v2","v2_m2")
resuid=as.data.frame(resuid)
library(GGally)
ggpairs(resuid)
ggpairs(part4[,-1])
# use fitted value to estimate mutant allele fraction
tmp.newAF=c(lm(mutect2~mutect,data=part1)$fitted.values,
      lm(mutect~varscan2,data=part2)$fitted.values,
      lm(mutect2~varscan2,data=part3)$fitted.values,
      lm(mutect~varscan2,data=part4)$fitted.values
      )
newAF = data.frame(index=rbind(part1,part2,part3,part4)[,1],TUMOR.AF.est=tmp.newAF)

D.norep.noback.coding_UTR.shared.wgs.withAF = merge(D.norep.noback.coding_UTR.wgs,newAF)
save(D.norep.noback.coding_UTR.shared.wgs.withAF, file="npc_passQC_wgs_updated_SYSUCC_final_20161030_object.RData")

tmp.af=withAF%>%filter(Gene_refGene%in%c("CYLD","TRAF3","TP53","NLRC5","SPATA3"))%>%select(SAMPLE,Gene_refGene,TUMOR.AF.est)
ggplot(tmp.af,aes(x=SAMPLE,y=Gene_refGene,fill=TUMOR.AF.est))+geom_tile()+
  scale_fill_gradient2()


### plot oncomap.
load("npc_passQC_wes_final_20161029_object.RData")
df = D.norep.noback.coding_UTR_intron.wes_wgs_HKSG.combined %>%
  filter(Func_refGene %in% c("exonic","splicing")) %>%
  select(sample=SAMPLE, gene=Gene_refGene, mutation=ExonicFunc_refGene,source=source)

# import mutsigCV results
mutsig = read.table("wes_exonic_intronic_20161030/wes_exonic_intronic_20161030.sig_genes.txt",
                    header=TRUE,sep="\t",stringsAsFactors = F)
mutsig$p[mutsig$p == 0] = 1e-16
mutsig$q[mutsig$q == 0] = 1e-16
# -log(P)
mutsig$log_P = round(-log(mutsig$p,base = 10),2)
mutsig$log_q = round(-log(mutsig$q,base = 10),2)

mutsig_plot = ggplot(head(mutsig,30), aes(x=gene, y=log_q)) +
  geom_bar(fill="skyblue", stat = "identity") + coord_flip() +
  scale_y_reverse(expand = c(0,0)) + geom_hline(yintercept = 3, lty=2) +
  theme(axis.line.y=element_blank()) 


jpeg("oncoplot.jpeg",width=1400,height=2000)
HM <- oncoplot(data=df,gene.annotation.plot = mutsig_plot,
         included.gene.list = c("TP53","CYLD","TRAF3","NFKBIA","NLRC5"),
         is.drop.gene=FALSE)
dev.off()

## only several genes
top30 = as.character(head(mutsig$gene,30))
oncoplot(data=df,gene.annotation.plot = mutsig_plot,
         included.gene.list = top30,
         is.drop.gene=FALSE,is.sort.gene = FALSE)

annot.table = unique(df %>% select(sample,source))

interested = c("TP53","CYLD","TRAF3","NFKBIA","NLRC5","RPL22","PRH2")
HM <- oncoplot(data=df, gene.annotation.plot = mutsig_plot,
               sample.annotation.table = annot.table,
               included.gene.list = interested,
               is.drop.gene=FALSE,is.sort.gene = FALSE)

# testing  start
test.df = D.norep.noback.coding_UTR_intron.wes_wgs_HKSG.combined %>%
  filter(Gene_refGene %in% interested) %>%
  select(sample=SAMPLE, gene=Gene_refGene, mutation=ExonicFunc_refGene)

sample_list=unique(test.df$sample)
test.annot.df = cbind(sample_list,
                      sample(c("T1","T2","T3"),length(sample_list),replace = T),
                      sample(c("WES","WGS","HK","SG"),length(sample_list),replace = T)
                      )
colnames(test.annot.df) = c("sample","TMN","source")
source("~/Documents/workspace/Github_sync_reposority/omics_scripts/oncoplot/oncoplot_heatmap.R")
oncoplot(data=test.df,gene.annotation.plot = mutsig_plot,
         sample.annotation.table = test.annot.df,
         is.drop.gene=FALSE,is.sort.gene = FALSE,is.sort.sample=FALSE)
### testing over

ggplot(onco, aes(x=sample,y=gene,fill=mutation)) + 
  geom_tile(width=1,height=1,colour="grey70") + 
  scale_fill_gradientn(colours=rainbow(4),na.value = "white") + 
  geom_text(aes(label=round(mutation,2)),angle=90) +
  xlab("") + ylab("") +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,size=10),axis.text.y=element_text(size=11))
  


### loading expression data
load("../../3-RNA-seq/expression_Hanbw/RNAseq_21NPCtumors_10inflam_Deseq_hanbw_cp_20161129.RData")
id_with_WTS=c("WES01080","WES01108","WES01172","WES01245","WES01293","WES01377","WES01393","WES01414",
              "WES01467","WES01482","WES01533","WGS01540","WGS01550","WGS01551","WGS01649","WES01652",
              "WES01655","WGS01675","WES01681","WES01694","WES01716","WES01641","WES01476","WES01111",
              "WES00952","WES01020")
names(id_with_WTS)=c("1080", "1108", "1172", "1245", "1293", "1377", "1393", "1414", 
                     "1467", "1482", "1533", "1540", "1550", "1551", "1649", "1652", 
                     "1655", "1675", "1681", "1694", "1662", "1641", "1476", "1111",
                     "952" , "1020")
 interested=c("TP53", "CYLD", "TRAF3", "NFKBIA", "NLRC5", "RPL22", "PRH2","TGFBR2","TET2")
tmp.c = colnames(RNAseq.expr.deseq.tumor)
tmp.c = gsub("WTS_(\\d+)_UN","\\1",tmp.c)
colnames(RNAseq.expr.deseq.tumor) = id_with_WTS[tmp.c]
tmp.expr = as.data.frame(cbind(RNAseq.expr.deseq.tumor[interested,],
                 RNAseq.expr.deseq.inflatmmation[interested,]
))
tmp.expr$gene = rownames(tmp.expr)
tmp.expr.df = melt(tmp.expr, variable.name = "sample", value.name = "expr")
# assign group
tmp.expr.df$group = "control"
tmp.expr.df$group[tmp.expr.df$sample %in% id_with_WTS] = "tumor"
tmp.expr.df$group = factor(tmp.expr.df$group,levels=c("tumor","control"))
pval=as.numeric(sapply(interested,function(x){
  try(with(tmp.expr.df %>% filter(gene==x),
           wilcox.test(expr~group))$p.value,silent = T)
}))
names(pval)=interested
tmp.expr.df$pval=sprintf("%.2e",pval[tmp.expr.df$gene])
ggplot(tmp.expr.df,aes(x=group,y=expr)) + 
  geom_boxplot(aes(fill=group)) + 
  geom_jitter() +
  geom_text(x=2,y=20000, aes(label=paste0("P=",pval)) ) +
  facet_wrap(~gene) +
  xlab("") + theme_bw() + theme(strip.text.x=element_text(size=12,face="bold"), 
                     axis.text.x=element_blank())

# mutation or not
tmp.mt = D.norep.noback.coding_UTR.shared.wes.withAF %>% 
  filter(Gene_refGene %in% interested) %>%
  filter(SAMPLE %in% id_with_WTS) %>%
  select(sample=SAMPLE,gene=Gene_refGene)
tmp.mt = as.data.frame(with(tmp.mt,table(sample,gene)))
tmp.mt$Freq[tmp.mt$Freq>1]=1
tmp.mt = tmp.mt %>%filter(Freq>0)
pval=as.numeric(sapply(interested,function(x){
  try(with(tmp.expr.df.tumor %>% filter(gene==x),
           wilcox.test(expr~mutation))$p.value,silent = T)
  }))
names(pval)=interested
tmp.expr.df.tumor = tmp.expr.df %>% filter(group=="tumor") %>%
  arrange(sample,gene)
tmp.expr.df.tumor$mutation = 0
for(i in 1:nrow(tmp.mt)){
  tmp.expr.df.tumor$mutation[tmp.expr.df.tumor$sample %in% tmp.mt$sample[i] &
                             tmp.expr.df.tumor$gene %in% tmp.mt$gene[i]] = 1
}
tmp.expr.df.tumor = within(tmp.expr.df.tumor,{
  mutation = factor(mutation,levels = c(0,1),labels = c("wildtype","mutation"))
})
tmp.expr.df.tumor$pval=round(pval[tmp.expr.df.tumor$gene],3)
ggplot(tmp.expr.df.tumor,aes(x=mutation,y=expr)) + 
  geom_boxplot(aes(fill=mutation)) + 
  geom_jitter() + 
  geom_text(x=2,y=20000, aes(label=paste0("P=",pval)) ) +
  facet_wrap(~gene) +
  xlab("") + theme_bw() + theme(strip.text.x=element_text(size=12,face="bold"), 
                                axis.text.x=element_blank())



# dn/ds plot plus gene expression
gene_summary = summarize_by_gene_annovar(D.norep.noback.coding_UTR.shared.wes.withAF %>%
                            filter(ExonicFunc_refGene != "non-exonic" & type=="SNV")  %>% 
                            filter(ExonicFunc_refGene != "unknown")
)
tmp.d=gene_summary%>%mutate(sample_frac=n_sample/474, dn_ds=(n_nonsense+n_missense)/n_silent)
tmp.max=max(tmp.d$dn_ds[!(tmp.d$dn_ds==Inf | is.nan(tmp.d$dn_ds))])
tmp.d$dn_ds[tmp.d$dn_ds==Inf]=tmp.max+5
tmp.d=tmp.d%>%filter(!is.nan(dn_ds)) %>%
  select(gene=Gene_refGene, sample_frac, dn_ds)

mean.RNAseq.expr.deseq.tumor$gene=rownames(mean.RNAseq.expr.deseq.tumor)
tmp.e = mean.RNAseq.expr.deseq.tumor %>% select(gene,median)

tmp.d = tmp.d %>% filter(gene %in% tmp.e$gene)
tmp.e = tmp.e %>% filter(gene %in% tmp.d$gene)
tmp.c = merge(tmp.d,tmp.e)



