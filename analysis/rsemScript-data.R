
rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"

ghpc <<- "aw30w@ghpcc06.umassrc.org"
lnc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.long_noncoding_RNAs.geneList"
pc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.annotation.pc.geneList"



RSEM_remap_doit <- function(){
  getDataTotalReadsBtwnReps_RSEM( reportFile=getFullPath("/data/rsemCapData-lpa-proc-REPORT.tab"),
                                  localDir = "/home/wespisea/data/RSEM-2/",
                                  remoteDir = file.path(rnaseqdir,"/rsem-hg19-gencodeV19"),
                                  outFile = getFullPath("/data/rsemCapData-v2-lpa-proc.tab"))  
}


RSEM_remap_doit_STAR <- function(){
  
  
  dummy <- getRSEMDataForOneCell( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                              localDir = "/home/wespisea/data/RSEM-fr/",
                                              remoteDir =  file.path(rnaseqdir,"/star-transcriptome/RSEM-fr"),
                                              getFile=TRUE,
                                              suffix=".genes.results")
  
  getDataTotalReadsBtwnReps_RSEM( reportFile=getFullPath("/data/rsemSTARmapped-lpa-proc-REPORT.tab"),
                                  localDir = "/home/wespisea/data/RSEM-fr/",
                                  remoteDir = file.path(rnaseqdir,"/star-transcriptome/RSEM-fr"),
                                  outFile = getFullPath("/data/rsemSTARmapped-v1-lpa-proc.tab"))  
  
  
  getReadDistro(localDir = "/home/wespisea/data/RSEM-fr-cnt/",
                            remoteDir = file.path(rnaseqdir,"/star-transcriptome/RSEM-fr"),
                            outFile = getFullPath("/data/rsemSTARmapped-v1-readDistro.tab"),
                            filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab")
}

RSEM_remap_doit_rsem <- function(){
  
  
  dummy <- getRSEMDataForOneCell( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                  localDir = "/home/wespisea/data/RSEM-self/",
                                  remoteDir =  file.path(rnaseqdir,"/rsem-hg19-gencodeV19"),
                                  getFile=TRUE,
                                  suffix=".genes.results")
  
  getDataTotalReadsBtwnReps_RSEM( reportFile=getFullPath("/data/rsemRsemmapped-lpa-proc-REPORT.tab"),
                                  localDir = "/home/wespisea/data/RSEM-self/",
                                  remoteDir = file.path(rnaseqdir,"/rsem-hg19-gencodeV19"),
                                  outFile = getFullPath("/data/rsemRsemmapped-v1-lpa-proc.tab"))  
  
  # re run
  getReadDistro( localDir = "/home/wespisea/data/RSEM-self-cnt/",
                 remoteDir = file.path(rnaseqdir,"rsem-hg19-gencodeV19"),
                 outFile = getFullPath("/data/rsemRsemmapped-v1-readDistro.tab"),
                filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab")
}

RSEM_remap_doit_bowtie <- function(){
  local.dir <- "/home/wespisea/data/RSEM-bowtie/"
  local.dir.reads <- "/home/wespisea/data/RSEM-bowtie-cnt/"
  
  remote.dir <- file.path(rnaseqdir,"/bowtie/RSEM")
  report.file <- getFullPath("/data/rsemBowtiemapped-lpa-proc-REPORT.tab")
  out.file <- getFullPath("/data/rsemBowtiemapped-v1-lpa-proc.tab")
  out.file.reads <- getFullPath("/data/rsemBowtiemapped-v1-readDistro.tab")
  
  dummy <- getRSEMDataForOneCell( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                  localDir = local.dir,
                                  remoteDir =  remote.dir,
                                  getFile=TRUE,
                                  suffix=".genes.results")
  
  getDataTotalReadsBtwnReps_RSEM( reportFile=report.file,
                                  localDir = local.dir,
                                  remoteDir = remote.dir,
                                  outFile =out.file)  
  
  
  getReadDistro( localDir =local.dir.reads,
                 remoteDir = remote.dir,
                 outFile = out.file.reads,
                 filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab")
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
getRSEMDataForOneCell <- function( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
                                   localDir = "/home/wespisea/data/flux/",
                                   remoteDir = file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14"),
                                   getFile=FALSE,
                                   suffix=".genes.results"){
 
  if(!file.exists(localDir)){
    dir.create(localDir)
  }
  
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & rnaExtract == "longPolyA")
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$remote <- paste0(remoteDir,"/",df.comb$bare,suffix)
  df.comb <- within(df.comb,{remoteFound <- hpc.file.exists(remote)})
  
  if(TRUE == getFile){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, "  ",localDir)
    write(o1,file="~/sandbox/rsemFetch")
    system("chmod u+x ~/sandbox/rsemFetch")
    suppressWarnings(system("~/sandbox/rsemFetch"))
    
  }
  df.comb$rsemFile <- paste0(localDir,df.comb$bare,".genes.results")
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"rsemFile", "bare")]
  colnames(df.comb) <- c("localization", "cell", "rnaExtract","replicate" ,"rsemFile", "bare")
  df.comb
}



readInRSEMGtfParsed <- function(file="/home/wespisea/data/RSEM/wgEncodeCshlLongRnaSeqHuvecNucleusPapFastqRep3.genes.results",spikeRet=FALSE){
  df <- read.csv(file=file,stringsAsFactors=FALSE ,sep="\t")
  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  #colnames(df) <- c("region", "gene_id", "reads", "length", "RPKM_byFluxC")
  df$region <- "other"
  df[which(df$gene_id %in% lnc), "region"] <- "lncRNA"
  df[which(df$gene_id %in% pc), "region"] <- "mRNA"
  df$isSpikeIn <- 0
  df[grep(df$gene_id, pattern = "ERCC"), "isSpikeIn" ] <- 1
  
  df
}



getTranscriptData_RSEM <- function(celltype,rnaExtract="longPolyA",cellMissing=c(),
                                   localDir = "/home/wespisea/data/flux/",
                                   remoteDir = file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14")){
  
  annot.df <- getRSEMDataForOneCell(filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",
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
  annot.df$fileExists <- sapply(annot.df$rsemFile,file.exists)
  annot.df$cellExtract = with(annot.df, paste0(cell,rnaExtract))
  counts <- ddply(annot.df, .(cellExtract), summarise, sum(fileExists))
  cellExtractFound <- counts[which(counts[,2] == 4),"cellExtract"]
  found.df <- annot.df[which(annot.df$cellExtract %in% cellExtractFound),]
  
  
  df.together <- data.frame()
  for(i in seq_along(found.df$rsemFile)){
    print(paste("finding data for -> ", found.df$cell[i]))
    df.local <- readInRSEMGtfParsed(file=found.df$rsemFile[i])
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

getDataTotalReadsBtwnReps_RSEM <- function(reportFile=getFullPath("/data/rsemCapData-lpa-proc-REPORT.tab"),
                                           localDir = "/home/wespisea/data/flux/",
                                           remoteDir = file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14"),
                                           outFile = getFullPath("/data/rsemCapData-lpa-proc.tab")){
  
  df.together <- getTranscriptData_RSEM( rnaExtract="longPolyA",
                                         localDir = localDir,
                                         remoteDir = remoteDir)
  print("got df.together") 
  report.df  <- as.data.frame(dplyr::group_by(df.together,cell,localization,replicate) %>%
                                dplyr::summarise(length(gene_id),
                                          mean(TPM),
                                          sum(TPM),
                                          mean(FPKM),
                                          sum(FPKM),
                                          sum(FPKM > 0),
                                          sum(expected_count)))
  report.df$experiment <- paste(ifelse(report.df$localization == "cytosol", "cyt", "nuc"),report.df$replicate,sep=".")
  colnames(report.df) <- c("cell", "localization", "replicate", "genesFound", "meanTPM", 
                           "sumTPM", "meanFPKM", "sumFPKM", "genesExpressed","sum_expected_count", "experiment")
  exportAsTable(df=report.df, file = reportFile)
  df.together <- as.data.frame(dplyr::group_by(df.together, cell, localization,rnaExtract,replicate) %>% 
                                 dplyr::mutate(FPKM_80norm = apply80norm(FPKM) * 1000000))
  
  #  group_by(df.together, cell, localization,rnaExtract,replicate) %>% summarise(mean(RPKM_80norm/transTotalRPKM, na.rm=TRUE))
  
  exportAsTable(file=paste0(outFile,".all"), df=df.together)
  df.together$gene_type <- df.together$region
  df.abbrev <- df.together[ c("region","replicate", "gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn", "FPKM_80norm", "TPM","FPKM")]
  
  df.rep.1 <- subset(df.abbrev, replicate == 1)
  df.rep.2 <- subset(df.abbrev, replicate == 2)
  
  df.cyt.rep1 <- subset(df.rep.1, localization == "cytosol")
  df.cyt.rep2 <- subset(df.rep.2, localization == "cytosol")
  idVars <-  c("gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn","replicate","region")
  idVarsNorep <- c("variable","gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn","region")
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
  
  idVars <-  c("gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn","replicate","region")
  idVarsNorep <- c("variable","gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn","region")
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
  
}


readInRSEMCntFile <- function(cntFile="~/sandbox/test.cnt"){
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



getReadDistro <- function(localDir = "/home/wespisea/data/RSEM-2-cnt/",
                          remoteDir = file.path(rnaseqdir,"/rsem-hg19-gencodeV19/"),
                          outFile = getFullPath("/data/rsemCapData-v2-readDistro.tab"),
                          filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab"){
  if(!file.exists(localDir)){
    dir.create(localDir)
  }
  
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & rnaExtract == "longPolyA")
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
    write(o1,file="~/sandbox/rsemFetch")
    system("chmod u+x ~/sandbox/rsemFetch")
    suppressWarnings(system("~/sandbox/rsemFetch 2>&1 > /dev/null"))  
  }
  
  
  df.comb$fileExists <- sapply(df.comb$local,file.exists)
  #annot.df$cellExtract = with(annot.df, paste0(cell,rnaExtract))
  #counts <- ddply(annot.df, .(cellExtract), summarise, sum(fileExists))
  #cellExtractFound <- counts[which(counts[,2] == 4),"cellExtract"]
  found.df <- df.comb[which(TRUE == df.comb$fileExists),]
  
  
  df.together <- data.frame()
  for(i in seq_along(found.df$local)){
    print(paste("finding data for -> ", df.comb$cell[i]))
    df.local <- readInRSEMCntFile(cntFile=found.df$local[i])
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






