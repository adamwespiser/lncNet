library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"
transGene <<- "/home/wespisea/work/research/researchProjects/coexpr/lncNET/data/gencode.v19.annotation.transGene.space"

ghpc <<- "aw30w@ghpcc06.umassrc.org"
lnc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.long_noncoding_RNAs.geneList"
pc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.annotation.pc.geneList"
apply80norm2 <- function(s){
  s.pos <- s[which(s > 0)]
  lowerp <- quantile(s.pos,0.1)
  upperp <- quantile(s.pos,0.9)
  s/sum(s[which(s > lowerp & s < upperp)])
}


apply80norm <- function(s){
  s.pos <- s[which(s > 0)]
  lowerp <- quantile(s.pos,0.2)
  upperp <- quantile(s.pos,0.8)
  s/sum(s[which(s > lowerp & s < upperp)])
}

eXpress_remap_doit_bowtie <- function(){
  rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  local.dir <- "/home/wespisea/data/eXpress_bowtie/"
  remote.dir <- file.path(rnaseqdir,"bowtie/eXpress")
  out.file <- getFullPath("/data/eXpress-bowtie-lpa-proc.tab")
  
  geteXpressDataForOneCell( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                            localDir = local.dir,
                            remoteDir = remote.dir,
                            getFile=TRUE,
                            suffix=".results.xprs",
                            remoteSuffix="results.xprs")
  
  
  getDataTotalReadsBtwnReps_eXpress( reportFile=getFullPath("/data/eXpress-nofr-lpa-proc-REPORT.tab"),
                                     localDir = local.dir,
                                     remoteDir = remote.dir,
                                     outFile = out.file)  
}


eXpress_remap_doit_2 <- function(){
  rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  
  geteXpressDataForOneCell( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                        localDir = "/home/wespisea/data/eXpress_nofr/",
                                        remoteDir = file.path(rnaseqdir,"star-transcriptome2/eXpress"),
                                        getFile=TRUE,
                                        suffix=".results.xprs",
                                        remoteSuffix="results.xprs")
  
  
  getDataTotalReadsBtwnReps_eXpress( reportFile=getFullPath("/data/eXpress-nofr-lpa-proc-REPORT.tab"),
                                     localDir = "/home/wespisea/data/eXpress_nofr/",
                                     remoteDir = file.path(rnaseqdir,"star-transcriptome2/eXpress"),
                                     outFile = getFullPath("/data/eXpress-nofr-lpa-proc.tab"))  
}




eXpress_remap_doit <- function(){
  getDataTotalReadsBtwnReps_eXpress( reportFile=getFullPath("/data/eXpress-lpa-proc-REPORT.tab"),
                                  localDir = "/home/wespisea/data/eXpress/",
                                  remoteDir = file.path(rnaseqdir,"star-transcriptome/eXpress"),
                                  outFile = getFullPath("/data/eXpress-lpa-proc.tab"))  
}

applyPseudoValByVar2 <- function(value,var){
  qFun <- function(x){
    x <- x[which(x > 0)]
    as.numeric(quantile(x,0.05)) 
  }
  pseudo <- tapply(value ,var, qFun)  
  cts <- as.numeric(sapply(levels(var), function(x)pseudo[[x]]))
  value + cts[var] 
}
geteXpressDataForOneCell <- function( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                   localDir = "/home/wespisea/data/eXpress/",
                                   remoteDir = file.path(rnaseqdir,"star-transcriptome/eXpress"),
                                   getFile=FALSE,
                                   suffix=".results.xprs",
                                   remoteSuffix="results.xprs"){
  
  if(!file.exists(localDir)){
    dir.create(localDir)
  }
  
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$remote <- paste0(remoteDir,"/",df.comb$bare,"/",remoteSuffix)
  df.comb$local <- paste0(localDir,"/",df.comb$bare,suffix)
  
  if(TRUE == getFile){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, "  ",df.comb$local)
    write(o1,file="~/sandbox/eXpressFetch")
    system("chmod u+x ~/sandbox/eXpressFetch")
    suppressWarnings(system("~/sandbox/eXpressFetch"))
    
  }
  #df.comb$File <- paste0(localDir,df.comb$bare,suffix)
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"local", "bare")]
  colnames(df.comb) <- c("localization", "cell", "rnaExtract","replicate" ,"local", "bare")
  df.comb
}



readIneXpressGtfParsed <- function(file="/home/wespisea/data/eXpress//wgEncodeCshlLongRnaSeqA549CytosolPapFastqRep3.results.xprs",spikeRet=FALSE){
  tdf    <- read.csv(file=file,stringsAsFactors=FALSE ,sep="\t")
  tg.df <- read.csv(file=transGene,stringsAsFactors=FALSE, sep=" ",header=FALSE)
  colnames(tg.df) <- c("gene_id","trans")
  
  m.df <- merge(x=tdf,y=tg.df,by.x="target_id",by.y="trans")
  m.df$geneFactor <- factor(m.df$gene_id)
  m.df$tpm = ifelse(m.df$solvable, m.df$tpm,0)
  m.df$fpkm = ifelse(m.df$solvable, m.df$fpkm,0)
  tt <- m.df[c("gene_id","tpm","fpkm","uniq_counts","eff_counts")]
  
  df  <- as.data.frame(  dplyr::group_by(tt,gene_id) %>%  
                                dplyr::summarize(a= sum(fpkm),
                                          b= sum(tpm),
                                          c=sum(uniq_counts),
                                          d=sum(eff_counts)))

  df$fpkm <- df$a
  df$a <- NULL
  df$tpm <- df$b
  df$b <- NULL
  df$uniq_counts <- df$c
  df$c <- NULL
  df$eff_counts <- df$d
  df$d <- NULL
  #gg <- group_by(tt, gene_id)
  #ggsum <- summarize(gg, sum(fpkm),sum(tpm))
  
  
  #df <- ddply(tt,.(gene_id),summarize,a=mean(tpm),b=mean(fpkm))
 # df <- genes.df 
  
  

  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  #colnames(df) <- c("region", "gene_id", "reads", "length", "RPKM_byFluxC")
  df$region <- "other"
  df[which(df$gene_id %in% lnc), "region"] <- "lncRNA"
  df[which(df$gene_id %in% pc), "region"] <- "mRNA"
  #df$isSpikeIn <- 0
  #df[grep(df$gene_id, pattern = "ERCC"), "isSpikeIn" ] <- 1
  
  df
}



getTranscriptData_eXpress <- function(celltype,rnaExtract,cellMissing=c("A549","IMR90"),
                                   localDir = "/home/wespisea/data/eXpress/",
                                   remoteDir = file.path(rnaseqdir,"star-transcriptome/eXpress")){
  
  annot.df <- geteXpressDataForOneCell(filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                    localDir = localDir,
                                    remoteDir = remoteDir,
                                    getFile=FALSE)
  
  if (!missing(celltype)){
    print("gathering all cell types")
    annot.df <- annot.df[which(annot.df$cell == celltype),]
  }
  if(!missing(rnaExtract)){
    annot.df <- annot.df[which(annot.df$rnaExtract == rnaExtract),]
  }
  # annot.df <- annot.df[-which(annot.df$cell %in% cellMissing),]
  annot.df$fileExists <- sapply(annot.df$local,file.exists)
  annot.df$cellExtract = with(annot.df, paste0(cell,rnaExtract))
  counts <- ddply(annot.df, .(cellExtract), summarise, sum(fileExists))
  cellExtractFound <- counts[which(counts[,2] == 4),"cellExtract"]
  found.df <- annot.df[which(annot.df$cellExtract %in% cellExtractFound),]
  
  
  df.together <- data.frame()
  for(i in seq_along(found.df$local)){
    print(paste("finding data for -> ", found.df$cell[i]))
    df.local <- readIneXpressGtfParsed(file=found.df$local[i])
    df.local$cell <- found.df$cell[i]
    df.local$localization <- found.df$localization[i]
    df.local$rnaExtract <- found.df$rnaExtract[i]
    df.local$replicate <- ifelse(found.df$replicate[i] > 2, found.df$replicate[i] -2, found.df$replicate[i])
    if (i == 1){
      df.together <- df.local
    } else{
      df.together <- rbind(df.together,df.local)
    }
  }
  df.together
}

getDataTotalReadsBtwnReps_eXpress <- function(reportFile=getFullPath("/data/eXpressCapData-nofr-lpa-proc-REPORT.tab"),
                                           localDir = "/home/wespisea/data/eXpress_nofr//",
                                           remoteDir = file.path(rnaseqdir,"star-transcriptome/eXpress2"),
                                           outFile = getFullPath("/data/eXpressCapData-nofr-lpa-proc.tab")){
  
  if(identical(globalenv(), environment())){
       variables <- ls()
reportFile=getFullPath("/data/eXpressCapData-nofr-lpa-proc-REPORT.tab")
    localDir = "/home/wespisea/data/eXpress_nofr//"
    remoteDir = file.path(rnaseqdir,"star-transcriptome2/eXpress")
    outFile = getFullPath("/data/eXpressCapData-nofr-lpa-proc.tab")
   }
  
  df.together <- getTranscriptData_eXpress( rnaExtract="longPolyA",
                                         localDir = localDir,
                                         remoteDir = remoteDir)
  df.together <- as.data.frame(dplyr::group_by(df.together, cell, localization,rnaExtract,replicate) %>% 
                                 dplyr::mutate(FPKM_80norm = apply80norm(fpkm) * 1000000))
  
  print("got df.together") 
  report.df  <- as.data.frame(dplyr::group_by(df.together,cell,localization,replicate) %>%
                                dplyr::summarise(length(gene_id),
                                          mean(tpm),
                                          sum(tpm),
                                          mean(fpkm),
                                          sum(fpkm),
                                          sum(fpkm > 0),
                                          sum(uniq_counts),
                                          sum(eff_counts)))
  report.df$experiment <- paste(ifelse(report.df$localization == "cytosol", "cyt", "nuc"),report.df$replicate,sep=".")
  colnames(report.df) <- c("cell", "localization", "replicate", "genesFound", "meanTPM", 
                           "sumTPM", "meanFPKM", "sumFPKM", "genesExpressed",
                           "uniqMappedReads","estimatedTotalReads", "experiment")
  report.df$fracUniq<- with(report.df, uniqMappedReads/estimatedTotalReads)
  exportAsTable(df=report.df, file = reportFile)
   
  #  group_by(df.together, cell, localization,rnaExtract,replicate) %>% summarise(mean(RPKM_80norm/transTotalRPKM, na.rm=TRUE))
  
  exportAsTable(file=paste0(outFile,".all"), df=df.together)
  
  df.together$gene_type <- df.together$region
  df.abbrev <- df.together[ c("region","replicate", "gene_id","gene_type", "localization","rnaExtract","cell", "FPKM_80norm", "tpm","fpkm")]
  
  df.rep.1 <- subset(df.abbrev, replicate == 1)
  df.rep.2 <- subset(df.abbrev, replicate == 2)
  
  df.cyt.rep1 <- subset(df.rep.1, localization == "cytosol")
  df.cyt.rep2 <- subset(df.rep.2, localization == "cytosol")
  idVars <-  c("gene_id","gene_type", "localization","rnaExtract","cell","replicate","region")
  idVarsNorep <- c("variable","gene_id","gene_type", "localization","rnaExtract","cell","region")
  df.cyt.rep1.melt <- melt(df.cyt.rep1, id.vars = idVars)
  df.cyt.rep2.melt <- melt(df.cyt.rep2, id.vars = idVars)
  df.cyt <- merge(df.cyt.rep1.melt, df.cyt.rep2.melt, by = idVarsNorep,suffixes=c(".rep1", ".rep2"))
  df.cyt$expr <- paste(df.cyt$localization,df.cyt$rnaExtract)
  df.cyt$value.rep1 <- ifelse(is.na(as.numeric(df.cyt$value.rep1)), 0, as.numeric(df.cyt$value.rep1))
  df.cyt$value.rep2 <-  ifelse(is.na(as.numeric(df.cyt$value.rep2)), 0, as.numeric(df.cyt$value.rep2))
  
  df.cyt$value.rep1.pseudo <- applyPseudoValByVar2(value= df.cyt$value.rep1, var=df.cyt$variable)
  df.cyt$value.rep2.pseudo <- applyPseudoValByVar2(value = df.cyt$value.rep2 , var=df.cyt$variable)
  
  df.cyt$rep1.frac <- df.cyt$value.rep1/(df.cyt$value.rep1 + df.cyt$value.rep2)
  df.cyt$rep1.frac.pseudo <- df.cyt$value.rep1.pseudo/(df.cyt$value.rep1.pseudo + df.cyt$value.rep2.pseudo)
  
  df.cyt$rep2.frac <- df.cyt$value.rep2/(df.cyt$value.rep1 + df.cyt$value.rep2)
  df.cyt$rep2.frac.pseudo <- df.cyt$value.rep2.pseudo/(df.cyt$value.rep1.pseudo + df.cyt$value.rep2.pseudo)
  
  df.cyt$rep.ratio <- df.cyt$value.rep1/( df.cyt$value.rep2)
  df.cyt$rep.ratio.pseudo <- df.cyt$value.rep1.pseudo/(df.cyt$value.rep2.pseudo)
  
  df.cyt$value.ave <- (df.cyt$value.rep1 + df.cyt$value.rep2)/2
  
  
  df.nuc.rep1 <- subset(df.rep.1, localization == "nucleus")
  df.nuc.rep2 <- subset(df.rep.2, localization == "nucleus")
  
  idVars <-  c("gene_id","gene_type", "localization","rnaExtract","cell", "replicate","region")
  idVarsNorep <- c("variable","gene_id","gene_type", "localization","rnaExtract","cell", "region")
  df.nuc.rep1.melt <- melt(df.nuc.rep1, id.vars=idVars)
  df.nuc.rep2.melt <- melt(df.nuc.rep2, id.vars = idVars)
  df.nuc <- merge(df.nuc.rep1.melt, df.nuc.rep2.melt, by = idVarsNorep,suffixes=c(".rep1", ".rep2"))
  df.nuc$expr <- paste(df.nuc$localization,df.nuc$rnaExtract)
  df.nuc$value.rep1 <- ifelse(is.na(as.numeric(df.nuc$value.rep1)), 0, as.numeric(df.nuc$value.rep1))
  df.nuc$value.rep2 <-  ifelse(is.na(as.numeric(df.nuc$value.rep2)), 0, as.numeric(df.nuc$value.rep2))
  
  df.nuc$value.rep1.pseudo <- applyPseudoValByVar2(value= df.nuc$value.rep1, var=df.nuc$variable)
  df.nuc$value.rep2.pseudo <- applyPseudoValByVar2(value = df.nuc$value.rep2 , var=df.nuc$variable)
  
  df.nuc$rep1.frac <- df.nuc$value.rep1/(df.nuc$value.rep1 + df.nuc$value.rep2)
  df.nuc$rep1.frac.pseudo <- df.nuc$value.rep1.pseudo/(df.nuc$value.rep1.pseudo + df.nuc$value.rep2.pseudo)
  
  df.nuc$rep2.frac <- df.nuc$value.rep2/(df.nuc$value.rep1 + df.nuc$value.rep2)
  df.nuc$rep2.frac.pseudo <- df.nuc$value.rep2.pseudo/(df.nuc$value.rep1.pseudo + df.nuc$value.rep2.pseudo)
  
  df.nuc$rep.ratio <- df.nuc$value.rep1/( df.nuc$value.rep2)
  df.nuc$rep.ratio.pseudo <- df.nuc$value.rep1.pseudo/(df.nuc$value.rep2.pseudo)
  
  df.nuc$value.ave <- (df.nuc$value.rep1 + df.nuc$value.rep2)/2
  
  
  df.cytNuc <- merge(df.cyt,df.nuc,by=c("gene_id","cell","variable"),suffixes=c(".cyt",".nuc"))
  df.cytNuc.rbind <- rbind(df.cyt,df.nuc)
  
  exportAsTable(file=paste0(outFile,".merge"), df=df.cytNuc)
  exportAsTable(file=outFile, df=df.cytNuc.rbind)
  
  if(identical(globalenv(), environment())){
    curr <- ls()
    rm(list=curr[!curr %in% variables])
  }
  
}


readIneXpressCntFile <- function(cntFile="~/sandbox/test.cnt"){
  df <- read.csv(file=cntFile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
  df.cnt <- df[4:(dim(df)[1]),]
  colnames(df.cnt) <- c("multiplicity","count")
  df.cnt$multiplicity <- as.numeric(df.cnt$multiplicity)
  statVec <- as(unlist(strsplit(df[1,1],split=" ")),"numeric")
  df.cnt$notMapped <- statVec[1]
  df.cnt$mapped <- statVec[2]
  df.cnt$totalRead <- statVec[4]
  df.cnt  
}



getReadDistro <- function(localDir = "/home/wespisea/data/eXpress-2-cnt/",
                          remoteDir = file.path(rnaseqdir,"/-hg19-gencodeV19/"),
                          outFile = getFullPath("/data/CapData-v2-readDistro.tab"),
                          filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab"){
  if(!file.exists(localDir)){
    dir.create(localDir)
  }
  
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$remote <- paste0(remoteDir,"/",df.comb$bare,".stat/",df.comb$bare,".cnt")
  df.comb$local <- paste0(localDir,df.comb$bare,".cnt")
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"remote","local", "bare")]
  colnames(df.comb) <- c("localization", "cell", "rnaExtract","replicate","remote","local", "bare")
  
  
  
  if(TRUE){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, "  ",localDir)
    write(o1,file="~/sandbox/Fetch")
    system("chmod u+x ~/sandbox/Fetch")
    suppressWarnings(system("~/sandbox/Fetch 2>&1 > /dev/null"))  
  }
  
  
  df.comb$fileExists <- sapply(df.comb$local,file.exists)
  #annot.df$cellExtract = with(annot.df, paste0(cell,rnaExtract))
  #counts <- ddply(annot.df, .(cellExtract), summarise, sum(fileExists))
  #cellExtractFound <- counts[which(counts[,2] == 4),"cellExtract"]
  found.df <- df.comb[which(TRUE == df.comb$fileExists),]
  
  
  df.together <- data.frame()
  for(i in seq_along(found.df$local)){
    print(paste("finding data for -> ", df.comb$cell[i]))
    df.local <- readIneXpressCntFile(cntFile=found.df$local[i])
    df.local$cell <- found.df$cell[i]
    df.local$localization <- found.df$localization[i]
    df.local$rnaExtract <- found.df$rnaExtract[i]
    df.local$replicate <- ifelse(found.df$replicate[i] > 2, found.df$replicate[i] -2, found.df$replicate[i])
    df.local$Exp<- paste0(df.local$cell,".",
                          ifelse(df.local$rnaExtract == "longNonPolyA","LNPA","lpa"),".",
                          df.local$replicate,".",
                          ifelse(df.local$localization=="cytosol","cyt","NUC" ) )
    
    df.together <- rbind(df.together,df.local)
    
  }
  
  exportAsTable(df=df.together, file = outFile)
  
}






