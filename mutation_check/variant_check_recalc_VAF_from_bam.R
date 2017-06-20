#######################################################################
# ToDo:  Check variants from BAM file, do hard-filtering and recalculate supporting
#           reads for alterate and reference alleles and VAF (variant allele frequency).
#        Assume the tumor and normal bam have a same naming style such as ${sample_id}_{T,B}.bam.
#
# Dependencies: Rsamtools, GenomicAlignments, stringr, Biostrings
#
# Input:  Argument 1: a mutation_file that should have the following columns with no header name.
#         But "#" line is allowed.
#           columns:
#           1: SAMPLE id
#           2: chromosome
#           3: position in 1-base coordination
#           4: reference allele
#           5: alterate allele
#           6: gene symbol (optional, not require)
#         Argument 2: output filenames
#
# Output: a table with the same lines and the first several columns as mutation_file, 
#         and add several columns after.
#           1: tumor.total.reads, total reads spanning the variant positions(including lowQual reads) in tumor.
#           2: tumor.pass.reads, number of reads that pass the read quality filters in tumor.
#           3: tumor.AD, reads supporting alterate alleles in tumor
#           4: tumor.RD, reads supporting reference alleles in tumor
#           5: tumor.VAF, variant allele frequency in tumor, VAF=AD/(AD+RD)
#           6: tumor.total.reads, total reads spanning the variant positions(including lowQual reads) in normal.
#           7: tumor.pass.reads, number of reads that pass the read quality filters in normal.
#           8: normal.AD, reads supporting alterate alleles in normal
#           9: normal.RD, reads supporting reference alleles in normal
#           10: normal.VAF, variant allele frequency in normal, VAF=AD/(AD+RD)
#           11: map.uniq, is the sequence around variant uniquely mapped to genome, thus less possibility of multiple
#               alignment or mapping error. TRUE or FALSE. ( Tue Jun 20 09:59:15 2017)
#           12: suppl.align, if map.uniq == FALSE, we could know how many regions are similar as the query one, 
#               and the details. ( Tue Jun 20 12:33:15 2017)
#
# Usage:  Rscript variant_check_recalc_VAF_from_BAM.R mutation_file output_file
#
# Author: Xiao-yu Zuo
# date: Thu Jun 15 10:03:25 2017
# 
# History:
# Tue Jun 20 09:59:15 2017 ------------------------------
#   Blast for sequences around variants.
#######################################################################
# Mon Jun 19 15:52:26 2017 ------------------------------
# suppress message
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(stringr))

# parse argument
args=commandArgs(TRUE)
mutation_file = as.character(args[1])
output_file = as.character(args[2])
path_to_bam = "/data/home/zuoxy/data/NPC/somatic/20160519_863closing/0-bam_GD/" # manual modify 
bam_surfix = ".bam" # manual modify

# set blast parameters and reference genome
ref_fasta = "/data/public/HomoSapiens/hg19/ucsc_hg19.fasta"
blast_dir = "/share/apps/blast/2.3.0/bin/" # path to blastn
blast_flank_bp = 150 # used to extract reference sequence of variants +- blast_flank_bp bps
blast_db = "/data/public/HomoSapiens/blastDB/hg19/ucsc_hg19" # blastdb path
nt = 2 # threads used in blastn

mutations = read.table(mutation_file, header=FALSE,sep="\t",stringsAsFactors = FALSE, comment.char = "#")
colnames(mutations) = c("sample","chr","pos","ref","alt","gene")[1:ncol(mutations)]

# Mon Jun 19 15:41:22 2017 ------------------------------
# in case of letting any column to be "logical"
col.logic = colnames(mutations)[sapply(mutations,class)=="logical"]
for(s in col.logic) mutations[,s] = str_sub(as.character(mutations[,s]),1,1)
rm(col.logic)


##############################
# Functions
##############################
#######################################################################
# Thu Jun 15 14:53:20 2017 ------------------------------
# check if a read is clipping by examing the cigar (S or H)
#######################################################################
isClipped <- function(cigar, tol=3){
  # if clipping occurs in both side, return TRUE for not confident
  # if clipping occurs in one side, but  > tol bases, return TRUE
  # if no clipping or only < tol bases clipping, return FALSE
  flag = FALSE
  if (str_detect(cigar,"[SH]")) {
    clip_str = str_extract_all(cigar,"[0-9]+[SH]",simplify = TRUE)
    if (ncol(clip_str) > 1) {# has both sided clipping
      flag = TRUE
    } else {
      clip_base = as.numeric(str_extract_all(clip_str,"[0-9]+", simplify=TRUE))
      if (clip_base > tol) { # one sided clipping but bases > tol
        flag = TRUE
      }
    }
  }
  return(flag)
}

isClipped_mateClipped <- function(cigar, rcigar, tol=3){
  # check for reads and their mates
  n1 = length(cigar)
  n2 = length(rcigar)
  if(n1 != n2) {
    stop("cigar and rcigar has unequal length! ")
  }
  
  isFail = sapply(1:n1, function(i){
    isClipped(cigar[i], tol) | isClipped(rcigar[i], tol)
  })
  
  return(isFail)
}

#######################################################################
#  check if a read is unimapped read, that is, no multiple alignment or
#  supplementary alignment found for this read. By examining the 'SA' or 
#  'XA' field
#######################################################################
isUniMapped <- function(SA){
  n = length(SA)
  f = !is.na(SA)
  if(length(f) == 1){
    return(rep(TRUE, n))
  }else{
    return(!f)
  }
}

#######################################################################
#  check if a read has sufficient mapping quality (MAPQ)
#######################################################################
isPass_MAPQ <- function(MAPQ, q=30){
  MAPQ >= q
}
isPass_MAPQ_mateMAPQ <- function(mq,rmq, q=30){
  isPass_MAPQ(mq, q) & isPass_MAPQ(rmq, q)
}

#######################################################################
# check if a read has correct insert size. (isize)
#######################################################################
isPass_isize <- function(isize, q=500){
  abs(isize) < q
}

#######################################################################
# check if a read carrying too many variants. (MD)
#######################################################################
isHyperMutatedRead <- function(MD, cigar, q=4){
  # Fri Jun 16 11:20:44 2017 ------------------------------
  # Some reads are weird as they have no "MD:Z" field in BAM.
  # Perhaps after indel realignmnet, something has changed?
  # for this case, I have to remove the reads with MD equals NA
  # and set the reads as non-pass.
  if(is.na(MD)){
    return(FALSE)
  }
  len.md = str_count(MD, "\\^*[A-Z]+")
  len.cigar = str_count(cigar,"I") # times of insertions
  len.nm = len.md + len.cigar
  return(len.nm >= q)
}

#######################################################################
# convert ASCII character to its Phred quality 
#######################################################################
convertPhredCharacter <- function(qual, phred = 33) {
  if(!phred %in% c(33, 64)) {
    stop("phred should be only 33 or 64 !")
  }
  as.numeric(charToRaw(qual)) - phred
}

#######################################################################
#  check if a read has (q<20) bases > n%, e.g. 50% (qual)
#######################################################################
isTooManyLowQualBases <- function(qual, phred = 33, q = 20, tol = 0.5) {
  # qual: bam QUAL field, the same as FASTQ QUAL field. ASCII character
  # phred: 33 or 64, for illumina Phred+33 or +64, respectively
  # q: threshold for base quality.
  # tol: maximun allowed threshold
  # return: TRUE for Yes, FALSE for No
  qual = convertPhredCharacter(qual, phred)
  isFail = (sum(qual < q)/length(qual)) >= tol
  return(isFail)
}

#######################################################################
#  check if a read carrying queried variants.
# Thu Jun 15 15:31:55 2017 ------------------------------
#   for this version, I have not yet check the exact sequence of insertion
# Fri Jun 16 13:09:08 2017 ------------------------------
#   add base quality (qual) filtering. If the mutated base is of low quality,
#   view it as not passed reads(i.e. let isMutated = FALSE)
#######################################################################
isMutatedRead <- function(start, mpos, ref, alt, MD, cigar ){
  # one read per input:
  # start: read start; 
  # mpos: mutation position; ref: ref base; alt: alt base
  # Note: pos, ref and alt should be 1-based coordination
  # MD: mismatch tag; cigar: cigar; 
  # return a list() with elements
  #   1. isMutated: TRUE: read carrying mutation, FALSE: read not carrying mutation
  #   2. mpos: the position of mutation on this read. 0 for notfound(non-mutated reads)
  type="SNP"
  if(ref == "-"){
    type = "INS"
  }else if(alt == "-"){
    type = "DEL"
  }
  if(is.na(MD)) {
    # Fri Jun 16 11:20:44 2017 ------------------------------
    # Some reads have no MD field (NA). See isHyperMutatedRead() for detail
    return(list(isMutated = FALSE, mpos = 0))
  }
  mp = mpos - start +1 # mutation point in a read
  len = end - start + 1 # length of a read
  clen = 0 # pos of mutation according for reference ( counts DEL)
  rlen = 0 # pos of mutation on a read (dont count DEL)
  isMutated = FALSE
  if(type %in% c("SNP","DEL")){
    mdz = str_match_all(MD,"[0-9]+\\^*[A-Z]*")[[1]]
    for(j in 1:nrow(mdz)){ # iterater for each mismatch
      x = mdz[j,]
      a = as.numeric(str_extract(x, "[0-9]+"))
      b = !is.na(str_extract(x, "\\^")) # is a DEL?
      c = str_extract(x,"[A-Z]+")
      if(is.na(c)){break} # the end of a read
      clen = clen + a
      rlen = rlen + a
      if(mp - clen == 1){ # meet the mutation
        if(type == "SNP" & !b & c == ref){
          isMutated = TRUE
          rlen = rlen + 1
          break
        }else if(type == "DEL" & b & c == ref){
          isMutated = TRUE
          rlen = rlen + 1
          break
        }
      }
      # if not meet the mutation, go on searching
      clen = clen + nchar(c)
      if(!b){# dont add rlen if DEL
        rlen = rlen + nchar(c)
      }
    }
  }else{
    # for INS, only count cigar
    # to update: need check the insertion bases equals to the altbases
    mdz = str_match_all(cigar,"[0-9]+[MID]")[[1]]
    for(j in 1:nrow(mdz)){ # iterater for each mismatch
      x = mdz[j,]
      a = as.numeric(str_extract(x, "[0-9]+"))
      c = str_extract(x,"[A-Z]+")
      if(c == "M"){
        clen = clen + a
        rlen = rlen + a
      }else if(c == "D"){
        clen = clen + a
      }else if(c == "I"){
        if(mp - clen == 1){
          if(a == nchar(alt)){
            isMutated = TRUE
            rlen = rlen + 1
            break
          }
        }
        rlen = rlen + a
      }
    }
  }
  return(list(isMutated = isMutated, mpos = ifelse(isMutated, rlen, 0)))
}

#######################################################################
# Thu Jun 15 14:01:12 2017 ------------------------------
# check and QC vriants from BAM file and recalulate the AD, RD and VAF
#######################################################################
scan_mutation_from_bam <- function(which, ref, alt, bamName, what, flag, tag, bq.cutoff = 20){
  # bq.cutoff: sequencing quality of nucleotide in FASTQ. bq < bq.cutoff would be removed
  AD = 0
  RD = 0
  VAF = 0
  TR = 0
  PR = 0  
  if(is.null(what)){
    flag = scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F,
                       isSecondaryAlignment = F, isNotPassingQualityControls = F,
                       isDuplicate = F)
  }
  if(is.null(flag)){
    what = c("qname","flag", "mapq", "isize", "seq", "qual")
  }
  if (is.null(tag)) {
    tag = c("MQ", "MC", "SA", "MD", "NM","XA")
  }
  
  # scan bam
  param <- ScanBamParam(flag=flag, which=which, what=what, tag=tag)
  bamCon = BamFile(bamName, asMates = FALSE) # see $mate_status for mated pairs (?BamFile)
  bam = as.data.frame(readGAlignments(bamCon,param=param,use.names = FALSE))
  
  TR = nrow(bam)
  if(TR < 1){
    return(c("total.read" = TR, "pass.read" = PR, "AD" = AD, "RD" = RD, "VAF" = VAF))
  }
  
  # filtering low confidence reads. TRUE for pass, FALSE for fail
  qc_flag = sapply(1:TR, function(i){
    flag = (!isClipped(bam$cigar[i], tol=5)) & 
      (isPass_MAPQ(bam$mapq[i])) &
      isPass_isize(bam$isize[i]) & 
      isUniMapped(bam$SA[i]) & 
      isUniMapped(bam$XA[i]) & 
      (!isHyperMutatedRead(bam$MD[i], bam$cigar[i], q = 4)) & 
      bam$width[i] > 60
    if(!is.na(bam$MC[i])) {
      base_flag = flag & (!isClipped(bam$MC[i], tol=8))
    }
    base_flag
  })
  
  bam.pass = bam[qc_flag & !is.na(qc_flag),]
  PR = nrow(bam.pass)
  
  # if no passed reads, return FALSE QC
  if(nrow(bam.pass) == 0) {
    return(c("total.read" = TR, "pass.read" = PR, "AD" = AD, "RD" = RD, "VAF" = VAF))
  } else {
    # is read carrying mutation?
    is.mutread = rep(F, nrow(bam.pass))
    mpos = rep(0, nrow(bam))
    for(j in 1:nrow(bam.pass)){
      z = isMutatedRead(bam.pass$start[j],pos,refbase,altbase,bam.pass$MD[j],bam.pass$cigar[j])
      is.mutread[j] = z[["isMutated"]]
      mpos[j] = z[["mpos"]]
    }
    # Fri Jun 16 13:23:47 2017 ------------------------------
    # if the mutated bases of read i have high base quality ? cannot check for DEL
    # if not pass, set is.mutread[i] = FALSE
    if(alt != "-") { # not check for DEL
      for(i in 1:length(is.mutread)) {
        if(!is.mutread[i]) { # skip non-mutated reads
          next
        }
        # extract qual characters. the same for SNP or INS
        # Fri Jun 16 15:46:16 2017 ------------------------------
        #   if cigar has H or S at the front, e.g. 2S134M, we should shift the 
        #   qual some bases (e.g. 2 for S, but 0 for H)
        to.shift = str_match(bam.pass$cigar[i], "[0-9]+[S]")[1]
        if(is.na(to.shift)) { # not found \\d+S, nothing change
          to.shift = 0
        }else{ # found S, get the shift bases
          to.shift = as.numeric(str_extract(to.shift, "[0-9]+"))
        }
        tmp.q = str_sub(bam.pass$qual[i], mpos[i] + to.shift, mpos[i] + to.shift + length(alt) - 1)
        tmp.q = convertPhredCharacter(tmp.q, phred = 33) # convert to number
        if(any(tmp.q < bq.cutoff)) { # assume Q20
          is.mutread[i] = FALSE
        }
      }
    }
    
    # recompute allele depth
    AD = sum(is.mutread)
    RD = sum(!is.mutread)
    VAF =AD/(AD + RD)
  }
  return(c("total.read" = TR, "pass.read" = PR, "AD" = AD, "RD" = RD, "VAF" = VAF))
}



#######################################################################
# Tue Jun 20 10:27:28 2017 ------------------------------
# BLAST to check mapping uniqueness of sequence.
#######################################################################
blast_check <- function(location, blast_dir, blast_db, blast_flank_bp, threads){
  # location: a data.frame to store position of variants, with colnames. 
  #           column 1: chr, column 2: pos
  # balst_dir: path to blastn
  # blast_db: makeblastdb prepared blast databases.
  # blast_flank_bp: how many basepairs around the variant should be extracted
  # threads: threads used to do blast
  # output: a nonduplicated data.frame with columns
  #         1. chr
  #         2. pos
  #         3. uniquely mapped or not(TRUE or FALSE)?
  suppressMessages(library(Biostrings))
  tmpfile = tempfile(pattern = "file", tmpdir = tempdir(), fileext = "") # use temp file
  
  location = unique(location) # make query uniqueness.
  # set ranges to extract sequences
  location$start = location$pos - blast_flank_bp
  location$end = location$pos + blast_flank_bp
  location$index = with(location, str_c(chr,":",pos))
  location$pos = NULL # remove pos column
  # build GRanges object
  gr = makeGRangesFromDataFrame(location, ignore.strand=TRUE, keep.extra.columns=TRUE)
  
  ## get fasta and write to file
  fa = FaFile(ref_fasta) # load indexed fasta
  seq = getSeq(fa, gr)
  names(seq) = gr$index # rename seqs
  # write tmp fasta file
  writeXStringSet(seq, tmpfile, format = "fasta")
  
  ## do blast
  cat("Doing blastn, please wait for a long? time....\n")
  include_flags = "qseqid qlen sacc slen qstart qend sstart send bitscore score pident nident mismatch gaps evalue qcovs"
  # set command line string
  cmd = sprintf("%s/blastn -query %s -db %s -outfmt '6 %s' -num_threads %i -max_target_seqs 5 -qcov_hsp_perc 50",
                blast_dir, tmpfile, blast_db, include_flags, nt)
  blast_result <- read.table(pipe(cmd)) # pipe
  # add column name
  colnames( blast_result ) <- as.vector(str_split(include_flags," ",simplify = T))
  cat("Blastn done! Now filtering...\n")
  
  ## check multiple alignments
  output = location
  output$map.uniq = FALSE
  output$suppl.align = "" # supplemental alignment
  # foreach location
  for(i in 1:nrow(output)){
    idx = output[i,"index"]
    sub.out = blast_result[blast_result$qseqid == idx, ] # get blast results subset
    nm = nrow(sub.out) # how many match records?
    if(nm == 0){ # nothing found
      stop(str_c(idx," have nothing match!!!!"))
    } else if (nm == 1){ # uniquely mapped
      if(sub.out$qend - sub.out$qstart == 2*blast_flank_bp) {
        output[i,"map.uniq"] = TRUE
      } else {
        print(sub.out)
        stop(str_c(idx," have unknown problem !!!"))
      }
    } else if (nm > 1) { # multiple alignment found
      # remove itself match
      sub.out = subset(sub.out, !(sstart == output[i,"start"] & send == output[i,"end"] & sacc == output[i,"chr"]))
      output[i, c("map.uniq","suppl.align")] = with(sub.out, {
        is.sa = ((qstart <= blast_flank_bp / 2 & qend >= blast_flank_bp * 1.5)  & # read length
            pident > 90) # similarity
        sa = ""
        if(any(is.sa)){
          sa = with(sub.out[is.sa,], str_c(
            sprintf("%s:%i-%i(%.1f)", sacc, sstart, send, pident), 
            collapse = ";"
          ))
        }
        c(!any(is.sa), sa)
      })
    }
  }
  rm(tmpfile)
  ## return 
  cat("Filtering done, return...\n")
  print(head(output))
  return(output)
}



#####################################################
mutations$tumor.total.reads = 0
mutations$tumor.pass.reads = 0
mutations$tumor.AD = 0
mutations$tumor.RD = 0
mutations$tumor.VAF = 0
mutations$normal.total.reads = 0
mutations$normal.pass.reads = 0
mutations$normal.AD = 0
mutations$normal.RD = 0
mutations$normal.VAF = 0

ib = 1 # chunk id
chunksize = 50 # chunk size to write output 

flag = scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F,
                   isSecondaryAlignment = F, isNotPassingQualityControls = F,
                   isDuplicate = F)
what = c("qname","flag", "mapq", "isize", "seq", "qual")
tag = c("MQ", "MC", "SA", "MD", "NM","XA")

cat("Scanning BAM files start !\n")

for(i in 1:nrow(mutations)){
  sample = mutations[i,1]
  chr = mutations[i,2]
  pos = mutations[i,3]
  start = pos - 0
  end = pos + 0
  refbase = mutations[i,4]
  altbase = mutations[i,5]
  which = GRanges(chr,IRanges(start, end))
  # for tumor checking
  bamName = str_c(path_to_bam, "/", sample, "_T", bam_surfix)
  tmp = scan_mutation_from_bam(which, refbase, altbase, bamName, what, flag, tag, bq.cutoff = 10)
  mutations[i, c("tumor.total.reads", "tumor.pass.reads", "tumor.AD", "tumor.RD", "tumor.VAF")] = tmp
  # for normal checking
  bamName = str_c(path_to_bam, "/", sample, "_B", bam_surfix)
  tmp = scan_mutation_from_bam(which, refbase, altbase, bamName, what, flag, tag, bq.cutoff = 10)
  mutations[i, c("normal.total.reads", "normal.pass.reads", "normal.AD", "normal.RD", "normal.VAF")] = tmp
 
  # logging
  if(i %% chunksize == 0){
    if(ib %% 10 == 0) {
      cat("\n") # each 10 "+" in a row
    }
    cat("+ ")
    ib = ib + 1
  }
  
}

cat("Scanning BAM files done!\nNow checking mapping uniqueness....\n")

## Tue Jun 20 01:56:11 2017 ------------------------------
## check uniqueness of sequence of +- flk bps around mutations,
## use BLAT here is also OK, but BLAST is used here.
#  Note that the returned data.frame is non-duplicated.
blast_check_result = blast_check(mutations[,c("chr","pos")],blast_dir,blast_db,blast_flank_bp,nt)

# combine 
mutations$index = with(mutations, str_c(chr, ":", pos))
rownames(blast_check_result) = blast_check_result$index
mutations[index, c("map.uniq", "suppl.align")] = blast_check_result[mutations$index, c("map.uniq", "suppl.align")]
mutations[,"index"] = NULL # remove index

cat("Checking mapping uniqueness done!\nNow writing output... \n")

write.table(mutations,
            file = output_file,
            sep = "\t",
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cat("all done! \n")



