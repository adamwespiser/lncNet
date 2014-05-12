
# getRpkmFromBamDataForOneCell(writeCopyScript=TRUE)
# $ chmod u+x ~/sandbox/rpkmFromBamFetch
# $ ~/sandbox/rpkmFromBamFetch
# processGeneToTranscript()

rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"

ghpc <<- "aw30w@ghpcc06.umassrc.org"
lnc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.long_noncoding_RNAs.geneList"
pc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.annotation.pc.geneList"


convertTransToGeneGtf <- function(transFile,geneFile){
  tf <- tempfile()
  system( paste("cat",transFile," | sed 's/[;\"]//g' |awk -F' ' '{print $10,$12,$14,$4,$5,$6}' > ",tf))
  trans.df <- read.csv(file=tf, sep=" ", stringsAsFactors=FALSE,header=FALSE)
  file.remove(tf)
  colnames(trans.df) <- c("gene_id", "transcript_id", "RPKM","startPos","stopPos","reads")
  trans.df$length <- with(trans.df, stopPos-startPos)
  trans.df$RPKM = as.numeric(trans.df$RPKM)
  trans.df$RPKM_tieBreaker <- trans.df$RPKM + runif(seq_along(trans.df$RPKM))/(10^9)
  gene.df <- as.data.frame(group_by(trans.df, gene_id) %.% filter(RPKM_tieBreaker == max(RPKM_tieBreaker))) 
  gene.df$RPKM_tieBreaker <- NULL
  gene.df$startPos <- NULL
  gene.df$stopPos <- NULL
  exportAsTable(df=gene.df  ,file=geneFile)
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
getRpkmFromBamDataForOneCell <- function( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",writeCopyScript=FALSE){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  #/project/umw_zhiping_weng/wespisea/rna-seq//starSpikeIn/
  df.comb$remote <- file.path(rnaseqdir,"starSpikeIn/",paste0(df.comb$bare,".trans.gtf"))
  
  if(writeCopyScript){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, " /home/wespisea/data/rpkmFromBam/",paste0(df.comb$bare,".trans.gtf"))
    write(o1,file="~/sandbox/rpkmFromBamFetch")
  }
  df.comb$rpkmFromBamFile <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".trans.gtf")
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"rpkmFromBamFile", "bare")]
  df.comb$rfbGene <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".genes.gtf")
  colnames(df.comb) <- c("localization", "cell", "rnaExtract", "replicate", 
                        "rpkmFromBamFile", "bare", "rfbGene")
  df.comb
}

processGeneToTranscript <- function(force=FALSE){
  annot.df <- getRpkmFromBamDataForOneCell()
  genes <- annot.df$rfbGene
  trans <- annot.df$rpkmFromBamFile
  for (i in seq_along(genes)){
    if((file.exists(trans[i]) && ! file.exists(genes[i])) || force == TRUE){
      print(paste(i,"/",length(genes),"file -> ",annot.df$bare[i]))
      convertTransToGeneGtf(transFile = trans[i], geneFile = genes[i])
    }
    
  }
  
  
}



readInrpkmFromBamGtfParsed <- function(file="/home/wespisea/data/rpkmFromBam/wgEncodeCshlLongRnaSeqHuvecNucleusPapFastqRep3.trans.gtf",
                                       spikeRet=FALSE){
  df <- read.csv(file=file,stringsAsFactors=FALSE ,sep="\t")
  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  #colnames(df) <- c("region", "gene_id", "reads", "length", "RPKM_byFluxC")
  df$region <- "other"
  df[which(df$gene_id %in% lnc), "region"] <- "lncRNA"
  df[which(df$gene_id %in% pc), "region"] <- "mRNA"
  df$isSpikeIn <- 0
  df[grep(df$gene_id, pattern = "ERCC"), "isSpikeIn" ] <- 1
  sumSpikeInRPKM <- sum(as.numeric(df[which(df$isSpikeIn == 1),"RPKM"]))/1000
  df$spikeIn_norm <- as.numeric(df$RPKM)/sumSpikeInRPKM
  df$readsPerKb <- with(df, (reads*10^3/length))
  
  spike.df <- getSpikeInDf()
  
  
  df.spike <- df[grep(pattern="ERCC",df$gene_id),]
  df.spike$readsPerKb <- df.spike$reads / (df.spike$length/1000)
  #df.spike$RPKM <- df.spike$readsPerKb / millionsOfReads
  spike <- merge(df.spike,spike.df, by="gene_id")
  #s.lm <- glm(Pool14nmol.ul ~ reads + length, data=spike,family="poisson")
  s.lm <- glm(Pool14nmol.ul ~ readsPerKb + 0, data=spike)
  df$concBySpikeIn <- predict(s.lm, newdata=df)
  
  df
}



getTranscriptData_rpkmFromBam <- function(celltype,rnaExtract){
  annot.df <- getRpkmFromBamDataForOneCell()
  if (!missing(celltype)){
    print("gathering all cell types")
    annot.df <- annot.df[which(annot.df$cell == celltype),]
  }
  if(!missing(rnaExtract)){
    annot.df <- annot.df[which(annot.df$rnaExtract == rnaExtract),]
  }
  #annot.df <- annot.df[-which(annot.df$cell %in% c("A549","IMR90")),]
  
  df.together <- data.frame()
  for(i in seq_along(annot.df$rfbGene)){
    print(paste("finding data for -> ", annot.df$cell[i]))
    df.local <- readInrpkmFromBamGtfParsed(file=annot.df$rfbGene[i])
    df.local$cell <- annot.df$cell[i]
    df.local$localization <- annot.df$localization[i]
    df.local$rnaExtract <- annot.df$rnaExtract[i]
    df.local$replicate <- ifelse(annot.df$replicate[i] > 2, annot.df$replicate[i] -2, annot.df$replicate[i])
    if (i == 1){
      df.together <- df.local
      
    } else{
      df.together <- rbind(df.together,df.local)
    }
  }
  df.together
  
}

getDataTotalReadsBtwnReps_rpkmFromBam <- function(){
  df.together <- getTranscriptData_rpkmFromBam( rnaExtract="longPolyA")
  df.together$RPKM <- as.numeric(df.together$RPKM)

  report.df  <- as.data.frame(group_by(df.together,cell,localization,replicate) %.%
                                summarise(length(gene_id),
                                          mean(RPKM),
                                          sum(RPKM),
                                          sum(RPKM > 0)))
  report.df$experiment <- paste(ifelse(report.df$localization == "cytosol", "cyt", "nuc"),report.df$replicate,sep=".")
  colnames(report.df) <- c("cell", "localization", "replicate", "genesFound", "meanRPKM", 
                           "sumRPKM","genesExpressed", "experiment")
  exportAsTable(df=report.df, file = getFullPath("/data/rpkmFromBAMCapData-lpa-proc-REPORT.tab"))
  
  df.together <- as.data.frame(group_by(df.together, cell, localization,rnaExtract,replicate) %.% 
                                 mutate(RPKM_80norm = apply80norm(RPKM) * 1000000))
  
  #  group_by(df.together, cell, localization,rnaExtract,replicate) %.% summarise(mean(RPKM_80norm/transTotalRPKM, na.rm=TRUE))
  
  exportAsTable(file=getFullPath("/data/rpkmFromBamDataAllCells.tab"), df=df.together)
  df.together$gene_type <- df.together$region
  df.abbrev <- df.together[ c("region","replicate", "gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn", "RPKM_80norm","RPKM","reads","concBySpikeIn","spikeIn_norm")]
  
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
  
  
  df.cytNuc <- rbind(df.cyt,df.nuc)
  #df.cytNuc[which(df.cytNuc$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc[which(df.cytNuc$gene_id %in% lnc),"region"] <- "lncRNA"
  
  
  
  exportAsTable(file=getFullPath("/data/rpkmFromBamCapData-lpa-proc.tab"), df=df.cytNuc)
  
}

