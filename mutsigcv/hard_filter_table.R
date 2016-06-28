###########################################################
# hard filter variants from tab delimited table
# this script will parse AD field, and get two
# new variable call DP and AD(only alt read depth)
###########################################################
args=(commandArgs(TRUE))
in_file=args[1]
out_file=args[2]
expr=args[3]

pop_freq_col=c("ExAC_ALL","X1000g2015aug_all","Kaviar_AF","HRC_AF")
data=read.table(in_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
# parse tumor and normal 'AD' field
tumor.adp=matrix(as.numeric(do.call(c,strsplit(data$TUMOR.AD,","))),nc=2,byrow=T)
normal.adp=matrix(as.numeric(do.call(c,strsplit(data$NORMAL.AD,","))),nc=2,byrow=T)
# add tumor ALT allele depth and total depth.
data=within(data,{TUMOR.AD=tumor.adp[,2];TUMOR.DP=tumor.adp[,1]+tumor.adp[,2];NORMAL.AD=normal.adp[,2];NORMAL.DP=normal.adp[,1]+normal.adp[,2];})
for(s in pop_freq_col){
	data[,s]=as.numeric(data[,s])
	data[is.na(data[,s]),s]=0
}
sub_data=subset(data, eval(parse(text=expr)))
write.table(sub_data, file=out_file, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


