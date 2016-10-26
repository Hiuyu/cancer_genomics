library(dplyr)
library(dtplyr)
library(reshape2)
library(ggplot2)
library(VennDiagram)
library(GenomicRanges)
setwd("C:/Users/LibJu/workspace/myData/NPC/cancer_genomics/caller_comparison")
source("0-helper_function.R")
message("start")
mutect2<-as.data.frame(data.table::fread("mutect2_pass_annotate.tsv"))
mutect<-as.data.frame(data.table::fread("mutect_pass_annotate.tsv"))
varscan.snv<-as.data.frame(data.table::fread("varscan2_pass_snv_annotate.tsv"))
varscan.indel<-as.data.frame(data.table::fread("varscan2_pass_indel_annotate.tsv"))
varscan=rbind(varscan.snv,varscan.indel)
## get the population frequency colnames
pop.freq.col=colnames(mutect2)[sapply(colnames(mutect2),function(x){
  grepl("ExAC",x)|grepl("1000g",x)
})]
pop.freq.col=c(pop.freq.col,"SYSUCC_CANCER_FREE_AF")
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
varscan$ExonicFunc_refGene=gsub("\\s+","_",varscan$ExonicFunc_refGene)
amb=with(varscan, wrapper_remove_ambigious_find_nearest(Gene_refGene,Func_refGene,ExonicFunc_refGene,GeneDetail_refGene,sep1=";",sep2="="))
varscan=within(varscan,{
  source=decide_source(SAMPLE)
  type=decide_variant_type(REF,ALT)
  Gene_refGene=amb$gene
  Func_refGene=amb$func_gene
  ExonicFunc_refGene=amb$func_exonic
  GeneDetail_refGene=amb$dist
})
save.image("raw_all_input_3caller_3callers.RData")
####### part 1 WES SNV comparison #####
D.mt2=mutect2%>%select(CHROM,POS,REF,ALT,SAMPLE,Gene_refGene,Func_refGene,ExonicFunc_refGene,TUMOR.GT,TUMOR.AD,NORMAL.GT,NORMAL.AD,one_of(pop.freq.col),SYSUCC_CANCER_FREE_NC,type,source,GeneDetail_refGene,AAChange_refGene)
D.mt=mutect%>%select(CHROM,POS,REF,ALT,SAMPLE,Gene_refGene,Func_refGene,ExonicFunc_refGene,TUMOR.GT,TUMOR.AD,NORMAL.GT,NORMAL.AD,one_of(pop.freq.col),SYSUCC_CANCER_FREE_NC,type,source,GeneDetail_refGene,AAChange_refGene)
D.vs=varscan%>%select(CHROM,POS=Start,REF,ALT,SAMPLE,Gene_refGene,Func_refGene,ExonicFunc_refGene,TUMOR.GT,TUMOR.AD,NORMAL.GT,NORMAL.AD,one_of(pop.freq.col),SYSUCC_CANCER_FREE_NC,type,source,GeneDetail_refGene,AAChange_refGene)
# total number of SNV count
print(paste("ALL variant count for mutect2, mutect, varscan is ",nrow(D.mt2), nrow(D.mt), nrow(D.vs)))
caller=c(rep("mutect2",nrow(D.mt2)),
         rep("mutect",nrow(D.mt)),
         rep("varscan2",nrow(D.vs))
)
D.all=cbind(rbind(D.mt2,D.mt,D.vs),caller)
num1=table(D.all$ExonicFunc_refGene,D.all$caller)
#ggplot(D.all%>%filter(ExonicFunc_refGene!="."),aes(x=ExonicFunc_refGene,fill=caller))+geom_bar(position="dodge")

## recode population frequency to be numeric
for(x in c(pop.freq.col,"SYSUCC_CANCER_FREE_NC")){
  D.all[D.all[,x]==".",x]="0"
  D.all[,x]=as.numeric(D.all[,x])
  D.all[is.na(D.all[,x]),x]=0
}
# very time consuming ! ! !
# find the max pop.freq in 1000G and ExAC
tmp.freq=D.all[,pop.freq.col]
D.all$max.pop.freq=0
print("finding max pop freq")
i=1
while(i<=nrow(D.all)){
  D.all$max.pop.freq[i]=max(tmp.freq[i,])
  if(i%%50000==0){message(paste0("parsed ",i/1000,"K sites"))}
  i=i+1
}
D.all$pop.freq=.assign_pop_freq_group(D.all$max.pop.freq)
D.all_snv_indel=D.all
# set index
D.all_snv_indel=D.all_snv_indel%>%mutate(index=paste0(SAMPLE,"_",CHROM,":",POS,"_",REF,">",ALT))
output_file="ALL_SNV_INDEL_3callers_20160926.RData"
save(D.all_snv_indel,file=output_file)
print(paste("all done, please load()",output_file,"in R"))

D.all_snv_indel$region="non-exon"
D.all_snv_indel$region[D.all_snv_indel$Func_refGene %in% c("exonic","UTR3","UTR5","splicing")] = "exon+UTR+splicing"
D.all_snv_indel=within(D.all_snv_indel,{ExonicFunc_refGene[ExonicFunc_refGene == "."]="non-exonic"})
tmp.df=with(D.all_snv_indel%>%filter(ExonicFunc_refGene!="."),prop.table(table(ExonicFunc_refGene,caller,region,type,source),2:5))
tmp.df=melt(tmp.df)
ggplot(tmp.df, aes(x=caller,y=value,fill=ExonicFunc_refGene))+geom_bar(stat="identity",position="stack")+ylab("frequency") +
  ggtitle("only Exonic variants") +facet_grid(type+region~source)

## generate concordance table for callers
tmp.rmcol=colnames(D.all_snv_indel)[grepl("AFR|EAS|EUR|AMR|FIN|OTH|SAS|NFE$",colnames(D.all_snv_indel),ignore.case = TRUE)]
D.all.1=D.all_snv_indel%>%select(-one_of(tmp.rmcol),-starts_with("TUMOR"),-starts_with("NORMAL"))
# reconstruct data
D.all.2=subset(D.all.1,!duplicated(D.all.1$index))
D.all.2$mutect=FALSE; D.all.2$mutect[D.all.2$index %in% D.all.1$index[D.all.1$caller=="mutect"]]=TRUE
D.all.2$mutect2=FALSE; D.all.2$mutect2[D.all.2$index %in% D.all.1$index[D.all.1$caller=="mutect2"]]=TRUE
D.all.2$varscan2=FALSE; D.all.2$varscan2[D.all.2$index %in% D.all.1$index[D.all.1$caller=="varscan2"]]=TRUE
D.all.2$share_times= D.all.2$mutect + D.all.2$mutect2 + D.all.2$varscan2

## filter germline
D.all.1.filter=D.all.1%>%filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3 )
D.all.2.filter=D.all.2%>%filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3 )
#D.all.2.exonic %>% filter(Func_refGene %in% c("UTR3","UTR5","splicing","exonic"))

## concordance analysis: venn plot
# all variants, before filtering
tmp.d=D.all.2
venn.diagram(
  list(
    mutect=which(with(tmp.d,mutect & type=="SNV")),
    mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
    varscan=which(with(tmp.d,varscan2 & type=="SNV"))
  ),
  main="All SNVs, before germline filtering",
  "venn_2caller_snv_all.tiff"
)
venn.diagram(
  list(
    mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
    varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
  ),
  main="All INDELs, before germline filtering",
  "venn_2caller_indel_all.tiff"
)
# all variants , after filtering
tmp.d=D.all.2
venn.diagram(
  list(
    mutect=which(with(tmp.d,mutect & type=="SNV")),
    mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
    varscan=which(with(tmp.d,varscan2 & type=="SNV"))
  ),
  main="All SNVs, after germline filtering",
  "venn_2caller_snv_filter.tiff"
)
venn.diagram(
  list(
    mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
    varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
  ),
  main="All INDELs, after germline filtering",
  "venn_2caller_indel_filter.tiff"
)
# exons only, before filtering
tmp.d=D.all.2%>%filter(region=="exon+UTR+splicing")
venn.diagram(
  list(
    mutect=which(with(tmp.d,mutect & type=="SNV")),
    mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
    varscan=which(with(tmp.d,varscan2 & type=="SNV"))
  ),
  main="All SNVs, before germline filtering",
  "venn_2caller_snv_all_exon.tiff"
)
venn.diagram(
  list(
    mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
    varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
  ),
  main="All INDELs, before germline filtering",
  "venn_2caller_indel_all_exon.tiff"
)
# exons only , after filtering
tmp.d=D.all.2.filter%>%filter(region=="exon+UTR+splicing")
venn.diagram(
  list(
    mutect=which(with(tmp.d,mutect & type=="SNV")),
    mutect2=which(with(tmp.d,mutect2 & type=="SNV")),
    varscan=which(with(tmp.d,varscan2 & type=="SNV"))
  ),
  main="All SNVs, after germline filtering",
  "venn_2caller_snv_filter_exon.tiff"
)
venn.diagram(
  list(
    mutect2=which(with(tmp.d,mutect2 & type=="INDEL")),
    varscan=which(with(tmp.d,varscan2 & type=="INDEL"))
  ),
  main="All INDELs, after germline filtering",
  "venn_2caller_indel_filter_exon.tiff"
)


save.image("workspace_call_comparison_20160926.RData")

tmp.t=D.all_snv_indel%>%filter(max.pop.freq > 0.05 & region=="non-exon")
tmp.t1=table(tmp.t$Gene_refGene,tmp.t$caller)


## repeat region (repeatmasker) overlap
load("/data/home/zuoxy/resource/genome/ucsc/hg19/annotation/Repeats/rmsk_GRange.RData")
a=with(D.all_snv_indel,convert_0base_to_1base(POS,REF,ALT))
D.all_snv_indel=within(D.all_snv_indel,{
  start=a$START
  end=a$END
  REF1=a$NEW.REF
  ALT1=a$NEW.ALT
  strand=rep("+",nrow(D.all_snv_indel))
})
D.all_snv_indel_grange = makeGRangesFromDataFrame(D.all_snv_indel,seqnames.field="CHROM",keep.extra.columns=T,starts.in.df.are.0based = TRUE)
tmp.rpid=subsetByOverlaps(D.all_snv_indel_grange,rmsk_grange,type="any",ignore.strand=T)
D.norep_snv_indel_grange=D.all_snv_indel_grange[!D.all_snv_indel_grange$index %in% tmp.rpid$index,]

## gene-level concordance
# raw data
tmp.d=D.all.2%>%filter(source=="WES" & Func_refGene!="intergenic")
tmp.d1=D.all.2 %>% filter(index %in% unique(D.norep_snv_indel_grange$index)) %>%
  filter(source=="WES" & Func_refGene!="intergenic") %>%
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_NC < 3)
venn.diagram(
  list(
    mutect=unique(tmp.d%>%filter(mutect)%>%select(Gene_refGene))[,1],
    mutect2=unique(tmp.d%>%filter(mutect2)%>%select(Gene_refGene))[,1],
    varscan2=unique(tmp.d%>%filter(varscan2)%>%select(Gene_refGene))[,1]
  ),
  main="All genes, after region & germline filtering",
  "venn_3caller_gene_filter.tiff"
)

## liuwsh heatmap gene comparison
liu.data=liu.data%>%mutate(index=paste0(SAMPLE,"_",Chr,":",Start,"-",End,"_",Ref,">",Alt))
D.all_snv_indel=as.data.frame(D.all_snv_indel_grange)
D.all_snv_indel=D.all_snv_indel%>%mutate(index=paste0(SAMPLE,"_",seqnames,":",start,"-",end,"_",REF1,">",ALT1))
tmp.d=D.all.2%>%filter(index%in%liu.data$index)
venn.diagram(
  list(
    mutect=unique(tmp.d%>%filter(mutect)%>%select(index))[,1],
    mutect2=unique(tmp.d%>%filter(mutect2)%>%select(index))[,1],
    varscan2=unique(tmp.d%>%filter(varscan2)%>%select(index))[,1]
  ),
  main="liu, variant",
  "liu.venn_3caller_variant.tiff"
)
tmp.d_grange=makeGRangesFromDataFrame(tmp.d,keep.extra.columns = T)
tmp.rp=subsetByOverlaps(tmp.d_grange,rmsk_grange,type="any",ignore.strand=T)


## top 50 high frequency gene.
z1.5=D.all.2%>%filter(share_times>0 & mutect2) %>%
  filter(source=="WES") %>% 
  filter(region!="non-exon") %>%
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>%
  group_by(Gene_refGene) %>% summarise(n=n()) %>% arrange(-n)

## only norep variants.
D.norep_snv_indel=as.data.frame(D.norep_snv_indel_grange)
D.norep.1=D.norep_snv_indel%>%select(-one_of(tmp.rmcol),-starts_with("TUMOR"),-starts_with("NORMAL"))
D.norep.2=subset(D.norep.1,!duplicated(D.norep.1$index))
D.norep.2$mutect=FALSE; D.norep.2$mutect[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="mutect"]]=TRUE
D.norep.2$mutect2=FALSE; D.norep.2$mutect2[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="mutect2"]]=TRUE
D.norep.2$varscan2=FALSE; D.norep.2$varscan2[D.norep.2$index %in% D.norep.1$index[D.norep.1$caller=="varscan2"]]=TRUE
D.norep.2$share_times= D.norep.2$mutect + D.norep.2$mutect2 + D.norep.2$varscan2
#save(D.norep.1,D.norep.2,file="norep_3callers_20160930.RData")

## exclude excluded or not-interested gene (PMID:22294350)
backlist_gene_PMID22294350=read.table("PMID22294350_Table_S7_gene_exclusion_list_final.txt",stringsAsFactors = F,header=F,sep="\t")
backlist_gene_PMID22294350=toupper(backlist_gene_PMID22294350[,1])
#save(backlist_gene_PMID22294350,file="backlist_gene_PMID22294350_TableS7.RData")
# excluded backlist genes and MUC family,
D.norep.noback.exon.2=D.norep.2%>%filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
  filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
  filter(! grepl("\\bOR\\d+",Gene_refGene)) %>%
  filter(region=="exon+UTR+splicing")

## make gene summary for several data
# z1: WES 3 callers; z2: WES 2 callers
z1.summary = summarize_by_gene_annovar(D.norep.noback.exon.2 %>% filter((share_times>2 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
                                         filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>% 
                                         filter(source=="WES") %>% filter(Func_refGene%in%c("exonic","splicing"))
                                       )
z2.summary = summarize_by_gene_annovar(D.norep.noback.exon.2 %>% filter((share_times>1 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
                                         filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>%
                                         filter(source=="WES") %>% filter(Func_refGene%in%c("exonic","splicing"))
                                       )
write.csv(z1.summary,"WES_3caller.csv")
write.csv(z2.summary,"WES_2caller.csv")

z1.1=z1.summary[,1:2];z1.1$id=1:nrow(z1.1)
z2.1=z2.summary[,1:2];z2.1$id=1:nrow(z2.1)
order=rep(c("3caller","2caller"),c(nrow(z1.summary),nrow(z2.summary)))
z12=cbind(rbind(z1.1,z2.1),order)
z12.1=dcast(z12,Gene_refGene~order,value.var = "n_sample")
colnames(z12.1)=c("gene","n.2caller","n.3caller")
z12.1=z12.1%>%arrange(gene)
z12.2=dcast(z12,Gene_refGene~order,value.var = "id")
colnames(z12.2)=c("gene","order.2caller","order.3caller")
z12.2=z12.2%>%arrange(gene)
z12.3=merge(z12.1,z12.2)
z12.3=z12.3%>%arrange(order.2caller)
write.csv(z12.3,"WES_2-3caller.csv")

tmp.z=z12.3%>%filter(order.3caller <=100 & order.2caller<=100 )
ggplot(tmp.z,aes(order.2caller,order.3caller))+geom_point()+geom_text_repel(aes(label=gene))+
  geom_hline(yintercept = c(50,100),lty=2)+geom_vline(xintercept = c(50,100),lty=2) +
  geom_smooth()

order.3caller <= 100 | 



library(org.Hs.eg.db)
tmp.z12=head(z12.3,1000)
gene <- z12.3$gene
geneid <- select(org.Hs.eg.db, keys=gene, keytype="SYMBOL",
                 columns="ENTREZID")



### combine all filtering steps
D.norep.noback.coding_UTR.wes=D.norep.2 %>%
  filter(!Gene_refGene%in%backlist_gene_PMID22294350) %>% 
  filter(!grepl("\\bMUC\\d+",Gene_refGene) ) %>%
  filter(!grepl("\\bOR\\d+",Gene_refGene)) %>%
  filter(region=="exon+UTR+splicing") %>%
  filter((share_times>1 & type=="SNV") | (share_times>1 & type=="INDEL")) %>% 
  filter(max.pop.freq < 0.01 & SYSUCC_CANCER_FREE_AF < 0.01) %>% 
  filter(source=="WES") %>%
  dplyr::select(Chr=seqnames,
                Start=start,
                End=end,
                Ref=REF,
                Alt=ALT,
                Func.refGene=Func_refGene,
                Gene.refGene=Gene_refGene,
                GeneDetail.refGene=GeneDetail_refGene,
                ExonicFunc.refGene=ExonicFunc_refGene,
                AAChange.refGene=AAChange_refGene,
                SAMPLE=SAMPLE,
                type=type
  ) %>%
  arrange(Chr,Start)

# to maftools
tmp.d=D.norep.noback.coding_UTR.wes%>%select(-type)
prefix="wes"
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
npc=read.maf(maf = maf_file, removeSilent = FALSE, useAll = FALSE)
tiff("bbb.tif",width=2000,height=1500)
oncoplot(maf = npc, top = 100,removeNonMutated=FALSE,fontSize = 12)
dev.off()
npc.oncodrive=oncodrive(maf=npc,minMut = 5,pvalMethod = "combined")
plotOncodrive(res = npc.oncodrive, fdrCutOff = 0.05, useFraction = FALSE)



## prepare for next step: signature and prognostics
D.coding_splicing.wes=D.norep.noback.coding_UTR.wes %>% 
  filter(Func.refGene %in% c("exonic", "splicing")) %>%
  filter(type=="SNV")
save(D.coding_splicing.wes,file="WES_coding-splicing_2016-10-18.RData")

## sample-gene mutation matrix
tmp.g=npc@gene.summary%>%filter(MutatedSamples>10)%>%dplyr::select(Hugo_Symbol)
tmp.g=tmp.g$Hugo_Symbol
tmp.d=D.coding_splicing.wes%>%filter(Gene.refGene%in%tmp.g)
gene.mut_mat=dcast(tmp.d,SAMPLE~Gene.refGene, value.var = "type", fun.aggregate = length)
gene.mut_mat[gene.mut_mat>1]=1
write.csv(gene.mut_mat,file = "sample_gene_matrix.csv")



# ka/ks plot
tmp.d=z2.summary%>%mutate(sample_frac=n_sample/474, ka_ks=n_missense/n_silent)
tmp.max=max(tmp.d$ka_ks[!(tmp.d$ka_ks==Inf | is.nan(tmp.d$ka_ks))])
tmp.d$ka_ks[tmp.d$ka_ks==Inf]=tmp.max+5
tmp.d=tmp.d%>%filter(!is.nan(ka_ks))
tmp.d1=tmp.d%>%filter(Gene_refGene%in%c("CYLD","TP53","RYR2","TRAF3","NFKBIA"))
intersted=c("CYLD","TP53","RYR2","TRAF3","NFKBIA","TNFAIP3","PTEN","KRAS","NRAS","ARID1A","BAP1","KMT2C","KMT2D","ATM","UBN1")
nfkb=c("TNF","RAN","TNFAIP3","CYLD","BCL10","RIPK2","PRKCA","TNFRSF1A","NFKB1","TRAF6","MALT1","ERC1","RELA","IKBKG","IKBKB","BIRC2","XPO1","CHUK","NOD2","UBE2D3","NFKBIA","ATM","BTRC")
top2=c("RHOA","DCN","EP300","TGFBR2","DGKG","PIK3C2G","PIK3CA","PLCZ1","CACNA1A","CACNA1C","DUSP6","FGF6","FGFR3","HRAS","NRAS","NTF3","TGFBR2","TP53","TRAF2","FGF23","TAOK1")
library(ggrepel)
ggplot(tmp.d,aes(x=sample_frac, y=log(ka_ks,base = 2)))+geom_point()+
  geom_vline(xintercept = c(0.01,0.03),lty=2)+geom_hline(yintercept = log(c(0.25,1,4),base = 2),lty=2)+
  geom_text_repel(data=tmp.d%>%filter(((ka_ks>4|ka_ks<0.25)&sample_frac>0.03)|sample_frac>0.05|Gene_refGene%in%intersted),aes(label=Gene_refGene))
ggplot(tmp.d,aes(y=n_missense, x=n_silent))+geom_point()+
  geom_vline(xintercept = c(0.01,0.03),lty=2)+geom_hline(yintercept = c(1,8),lty=2)+
  geom_text_repel(data=tmp.d%>%filter(ka_ks>8&sample_frac>0.03),aes(label=Gene_refGene))
ggplot(tmp.d,aes(x=sample_frac, y=log(ka_ks,base = 2)))+geom_point()+
  geom_vline(xintercept = c(0.01,0.03),lty=2)+geom_hline(yintercept = log(c(0.25,1,4),base = 2),lty=2)+
  geom_text_repel(data=tmp.d%>%filter(Gene_refGene%in%top2),aes(label=Gene_refGene))

