###  the priority of functions, suitable only for ANNOVAR output
priority_gene_func<-function(){
  pri=c("exonic","splicing","UTR5","UTR3","intronic","upstream","downstream","ncRNA_exonic","ncRNA_splicing","ncRNA_intronic","intergenic")
  xx=1:length(pri)
  names(xx)=pri
  xx
}

###  the priority of functions, suitable only for ANNOVAR output
priority_exonic_func<-function(){
  pri=c(
    "stopgain" = 1,
    "stoploss" = 1,
    "frameshift_insertion" = 1,
    "frameshift_deletion" = 1,
    "nonframeshift_deletion" = 2,
    "nonframeshift_insertion" = 2,
    "nonsynonymous_SNV" = 2,
    "synonymous_SNV" = 3,
    "unknown" = 9,
    "." = 10
  )
  pri
}

### select most priority func and gene 
remove_ambigious<-function(gene,func_gene,func_exonic,sep=";"){
  gene.list=strsplit(gene,split = sep)
  func_gene.list=strsplit(func_gene,split = sep)
  func_exonic.list=strsplit(func_exonic,split = sep)
  pri_gene=priority_gene_func()
  pri_exonic=priority_exonic_func()
  score_gene.list=sapply(func_gene.list,function(x){
    score=pri_gene[x]
    return(which.min(score))
  })
  score_exonic.list=sapply(func_exonic.list,function(x){
    score=pri_exonic[x]
    return(which.min(score))
  })
  new.func_gene=sapply(1:length(func_gene.list),function(i){func_gene.list[[i]][score_gene.list[[i]]]})
  new.func_exonic=sapply(1:length(func_exonic.list),function(i){func_exonic.list[[i]][score_exonic.list[[i]]]})
  new.gene=sapply(1:length(func_gene.list),function(i){gene.list[[i]][score_gene.list[[i]]]})
  list(gene=new.gene,func_gene=new.func_gene,func_exonic=new.func_exonic)
}

### find gene nearest if and only if the func is intergenic
.find_nearest_gene_single<-function(gene, dist, sep1=";", sep2="="){
  if(!grepl("dist",dist)) stop("zzzzzzzz")
  gene.list=unlist(strsplit(gene,split = ","))
  dist.list=unlist(strsplit(dist,split = sep1))
  dist.list=gsub(paste0("dist",sep2),"",dist.list)
  id=NA
  if(all(dist.list=="NONE")){
    return(list(gene=NA,dist=NA))
  }else if(any(dist.list=="NONE")){
    id=which(dist.list!="NONE")
  }else{
    id=which.min(as.numeric(dist.list))
  }
  list(gene=gene.list[id],dist=as.numeric(dist.list[id]))
}

find_nearest_gene<-function(gene, dist,sep1=";",sep2="="){
  gene.dist.mat=sapply(1:length(gene),function(i){
    if(grepl("^dist",dist[i])){
      dd=.find_nearest_gene_single(gene[i],dist[i],sep1, sep2)
      return(c(dd$gene,dd$dist))
    }else{
      return(c(gene[i],dist[i]))
    }
  })
  list(gene=gene.dist.mat[1,],dist=gene.dist.mat[2,])
}

### parsing the GATK AD field to two columns: AD(ALT DEPTH) and DP(total Depth)
# require the input AD be character vector
# e.g. 25,2 -> 2  27
# e.g. 25,2,0,3 -> 5  30
# 2016-09-24
parse_GATK_AD_field<-function(AD){
  if(typeof(AD)=="character"){
    # do nothing
  }else if(typeof(AD)=="list"){
    AD=as.chracter(AD)
  }
  AD[is.na(AD)]="0,0" # recode NA to 0,0, means not covered
  n=length(AD)
  ad=numeric(n)
  dp=numeric(n)
  af=numeric(n)
  i=1 # looping over all sites
  while(i <= n){
    s=as.numeric(unlist(strsplit(AD[i],split = ",")))
    ad[i]=sum(s[-1]) ## sum over all read depth of non-ref alleles (Note: may not be precise for MNP and INDEL)
    dp[i]=sum(s)
    af[i]=ifelse(dp[i]==0,0,ad[i]/dp[i])
    if(i%%50000==0){message(paste0("parsed ",i/1000,"K sites"))}
    i=i+1
  }
  return(list(AD=ad,DP=dp,AF=af))
}

### decide the type of variants: SNV or INDEL
# 2016-09-24
decide_variant_type<-function(REF,ALT){
  if(typeof(REF)!="character" | typeof(ALT)!="character"){
    stop("the type of REF or ALT not correct. Only Character Vector are required")
  }else if(length(REF) != length(ALT)){
    stop("the length of REF not equal to ALT")
  }
  n=length(REF)
  type=rep("SNV",n)
  idx=(nchar(REF)!= nchar(ALT)) | REF=="-" | ALT=="-"
  type[idx]="INDEL"
  return(type)
}

### decide the source of variants: WES or WGS
# 2016-09-24
decide_source<-function(SAMPLE){
  if(typeof(SAMPLE)!="character"){
    stop("the type of SAMPLE not correct. Only Character Vector are required")
  }
  source=rep("UNKNOWN",length(SAMPLE))
  source[grepl("WES|wes",SAMPLE)]="WES"
  source[grepl("WGS|wgs",SAMPLE)]="WGS"
  return(source)
}

### wrapper for remove_ambigious and pickup nearest gene
# also see remove_ambigious() and find_nearest_gene()
# 2016-09-23
wrapper_remove_ambigious_find_nearest<-function(gene, func_gene, func_exonic, gene_detail, sep1, sep2){
  ## select most 'powerful' annotation for multi-hit gene
  amb=remove_ambigious(gene,func_gene,func_exonic,sep1)
  ## select nearest gene for intergenic
  amb1=find_nearest_gene(gene,gene_detail,sep1, sep2)
  amb$gene=amb1$gene
  amb$dist=amb1$dist
  return(amb)
}

### decide population frequency group
.assign_pop_freq_group<-function(freq){
  if(class(freq)!="numeric"){
    stop("FREQ require numeric format!")
  }
  freq_group=cut(freq,breaks = c(0,0.0005,0.001,0.005,0.01,0.05,0.1,1),
                 labels = c("<0.05%","0.05%~0.1%","0.1%~0.5%","0.5%~1%","1%~5%","5%~10%",">=10%"),
                 include.lowest = T, right=F
                 )
  freq_group=freq_group
  return(freq_group)
}




### 0-based coordination to 1-based , e.g. format from VCF -> MAF
# 2016-08-31
convert_0base_to_1base<-function(POS,REF,ALT){
  # POS: VCF 0-based  position
  # REF: VCF REF alleles
  # ALT: VCF ALT alleles
  ind.indel=(nchar(REF)>1 | nchar(ALT) >1)
  len.ref=nchar(REF)
  len.alt=nchar(ALT)
  new.REF=REF
  new.ALT=ALT
  START=POS
  new.REF[len.ref==1 & ind.indel]="-"
  new.ALT[len.ref==1 & ind.indel]=gsub("^\\w","",new.ALT[len.ref==1 & ind.indel])
  new.ALT[len.alt==1 & ind.indel]="-"
  new.REF[len.alt==1 & ind.indel]=gsub("^\\w","",new.REF[len.alt==1 & ind.indel])
  START[ind.indel]=START[ind.indel]+1
  END=START+(len.ref-1)
  END[len.ref>1]=END[len.ref>1]-1
  return(data.frame(START=START,END=END,NEW.REF=new.REF,NEW.ALT=new.ALT,stringsAsFactors = F))
}


#####  get the context and mutation type from reference genome ######
# 2016-09-30

######  decide mutation type ########
# 2019-09-30
### include 96 tri-context types and 6 1-based type, return a list()
.load_lego_type<-function(){
  # [ref][alt][5'context][3'context] or 
  # complement([ref][alt])reverseComplement([5'context][3'context])
  lego_96_type=c("TGTT", "TGGT", "TGCT", "TGAT", "TGTG", "TGGG", "TGCG", "TGAG",
                 "TGTC", "TGGC", "TGCC", "TGAC", "TGTA", "TGGA", "TGCA", "TGAA",
                 "TCTT", "TCGT", "TCCT", "TCAT", "TCTG", "TCGG", "TCCG", "TCAG",
                 "TCTC", "TCGC", "TCCC", "TCAC", "TCTA", "TCGA", "TCCA", "TCAA",
                 "TATT", "TAGT", "TACT", "TAAT", "TATG", "TAGG", "TACG", "TAAG",
                 "TATC", "TAGC", "TACC", "TAAC", "TATA", "TAGA", "TACA", "TAAA",
                 "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT",
                 "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT",
                 "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT",
                 "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT",
                 "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT",
                 "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT"
  )
  lego_6_type=rep(c("T>G","T>C","T>A","C>A","C>G","C>T"),2)
  names(lego_6_type)=c(
    c("T>G","T>C","T>A","C>A","C>G","C>T"), # the straight type
    c("A>C","A>G","A>T","G>T","G>C","G>A") # the complement type
  )
  list(lego_96_type=lego_96_type, lego_6_type=lego_6_type)
}

#####  define single-base mutation type ######
# 2016-09-30
##  [ref]>[alt] or complement([ref])>complement([alt])
## if input is a INDEL, return NA, and warning
decide_1base_mutation_type<-function(REF,ALT){
  require(Biostrings)
  # 6 types of mutation
  lego_6_type=.load_lego_type()[["lego_6_type"]]
  type=decide_variant_type(REF,ALT)
  is_SNV=type=="SNV"
  if(any(!is_SNV)){
    message(paste0("found ",sum(!is_SNV)," records not SNV! May be INDELs? ,e.g. "
                   ,paste0(head(which(!is_SNV)),collapse = ", "), 
                   ". return NA for these records"))
  }
  idx=paste0(REF[is_SNV],">",ALT[is_SNV])
  mutation_type=rep(NA,length(REF))
  mutation_type[is_SNV]=lego_6_type[idx]
  return(mutation_type)
}
















######  bayesNMF for mutation signature discovery ###########
### make lego matrix for BayesNMF









###### make gene mutation summarys #######
# 2016-09-30
summarize_by_gene_annovar<-function(D, add=NULL){
  # require data.frame, 
  # requireif colname: REF, ALT, Gene_refGene, Func_refGene, ExonicFunc_refGene, SAMPLE
  # the REF and ALT should be 1-based coding
  require(dplyr)
  require(reshape2)
  rq.col=c("REF", "ALT", "Gene_refGene", "Func_refGene", "ExonicFunc_refGene", "SAMPLE")
  rq.in=rq.col%in%colnames(D)
  # if any col not exists, stop
  if(any(!rq.in)){ 
    stop(paste0(
      paste0(rq.col[!rq.in],collapse = ", "),
      " not exists!")
      )
  }
  D=D[,rq.col]
  ## remove intergenic, or ncRNA variants
  D=D%>%filter(! Func_refGene %in% c("intergenic","ncRNA_exonic","nc_RNA_intronic", "ncRNA_splicing",
                                          "downstream","upstream")
  )
  D$type=decide_variant_type(D$REF,D$ALT)
  D$mutation_type=decide_1base_mutation_type(D$REF,D$ALT)
  # summarize base statistics
  summary_1=D%>%group_by(Gene_refGene)%>%
    summarise(n_sample=n_distinct(SAMPLE), 
              n_mutation=n(), 
              n_SNV=sum(type=="SNV"),
              n_indel=sum(type=="INDEL"),
              n_coding_SNV=sum(Func_refGene=="exonic" & type=="SNV"),
              n_UTR_splice_SNV=sum(Func_refGene %in% c("UTR3","UTR5","splicing") & type=="SNV"),
              n_intron_SNV=sum(Func_refGene=="intronic" & type=="SNV"),
              n_coding_indel=sum(Func_refGene=="exonic" & type=="indel"),
              n_UTR_splice_indel=sum(Func_refGene %in% c("UTR3","UTR5","splicing") & type=="indel"),
              n_intron_indel=sum(Func_refGene=="intronic" & type=="indel"),
              n_missense=sum(ExonicFunc_refGene=="nonsynonymous_SNV"),
              n_silent=sum(ExonicFunc_refGene=="synonymous_SNV"),
              n_nonsense=sum(ExonicFunc_refGene %in% c("stopgain","stoploss")),
              n_frameshift=sum(ExonicFunc_refGene %in% c("frameshift_deletion","frameshift_insertion")),
              n_nonframeshift=sum(ExonicFunc_refGene %in% c("nonframeshift_deletion","nonframeshift_insertion")),
              n_unknown=sum(ExonicFunc_refGene %in% c("unknown"))
              )
  # summarize mutation type
  tmp.s2=D%>%group_by(Gene_refGene,mutation_type) %>%
    summarise(value=n())
  tmp.s2$mutation_type[is.na(tmp.s2$mutation_type)]="INDEL"
  summary_2=dcast(tmp.s2,Gene_refGene~mutation_type)
  summary_2[is.na(summary_2)]=0 # assign 0 to NA
  if(any(colnames(summary_2)=="INDEL")){
    summary_2=summary_2%>%dplyr::select(-INDEL) # INDEL col removed
  }
  # combine different summary part
  summary_all=merge(summary_1,summary_2,by="Gene_refGene")
  summary_all=summary_all%>%arrange(-n_sample)
  return(summary_all)
}



