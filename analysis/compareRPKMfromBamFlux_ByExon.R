
getRpkmFromBamDataForOneCellByExon <- function( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",writeCopyScript=FALSE){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df$readLength <- as.numeric(gsub(sapply(strsplit(df$readType, "x"), function(x)x[2]),pattern="D",replacement=""))
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  #/project/umw_zhiping_weng/wespisea/rna-seq//starSpikeIn/
  df.comb$remote <- file.path(rnaseqdir,"starSpikeIn/",paste0(df.comb$bare,".transByExon.gtf"))
  
  if(writeCopyScript){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, " /home/wespisea/data/rpkmFromBam/",paste0(df.comb$bare,".transByExon.gtf"))
    write(o1,file="~/sandbox/rpkmFromBamExonFetch")
  }
  df.comb$rpkmFromBamFile <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".transByExon.gtf")
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"rpkmFromBamFile", "bare","read1.readLength")]
  df.comb$rfbTrans <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".transFromExon.gtf")
  colnames(df.comb) <- c("localization", "cell", "rnaExtract", "replicate", 
                         "rpkmFromBamFile", "bare", "readLength","rfbGene")
  df.comb
}


convertExonsToTrans <- function(transFile,exonFile){
  tf <- tempfile()
  system( paste("cat",exonFile," | sed 's/[;\"]//g' |awk -F' ' '{print  $4,$5,$6,$10,$12,$14}' > ",tf))
  exon.df <- read.csv(file=tf, sep=" ", stringsAsFactors=FALSE,header=FALSE)
  file.remove(tf)
  colnames(exon.df) <- c("startPos","stopPos","reads","gene_id", "transcript_id", "RPKM")
  exon.df$RPKM <- NULL
  exon.df$length <- with(exon.df,stopPos - startPos)
  exon.df$stopPos <- NULL
  exon.df$startPos <- NULL
  
  trans.df <- as.data.frame(group_by(exon.df, transcript_id) %.% 
                              summarise(reads = sum(reads),
                                        length = sum(length),
                                        gene_id = gene_id[1]))
  trans.df$readsPerLen <- with(trans.df, reads/length)
  
  exportAsTable(df=trans.df  ,file=transFile)
}

processCellsMaxTransExprReads <- function(){
  annot.df <- getRpkmFromBamDataForOneCellByExon()
  annot.df <- annot.df[which(annot.df$rnaExtract == "longPolyA"),]
  annot.df <- annot.df[-which(annot.df$cell == "H1-hESC"),]
  trans <- annot.df$rfbTrans
  exons <- annot.df$rpkmFromBamFile
  annot.df$rep <- ifelse(annot.df$replicate >2,annot.df$replicate -2,annot.df$replicate )
#  annot.df$transFullReads<- gsub(x=annot.df$rfbGene,pattern="genes",replacement="transFullReads")
 # transFullReads <- annot.df$transFullReads 
  sapply(seq_along(transFullReads), function(x)convertExonsToTrans(exonFile=exons[x],transFile=trans[x]))
  df.together <- data.frame()
  for ( cell in unique(annot.df$cell)){
    print(cell)
    a.cell <- annot.df[which(annot.df$cell == cell),]
    a.cell.cyt1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "cytosol"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt1$loc <- "cytosol"
    a.cell.cyt1$rep <- 1
    a.cell.nuc1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "nucleus"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc1$loc <- "nucleus"
    a.cell.nuc1$rep <- 1
    a.cell.cyt2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "cytosol"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt2$loc <- "cytosol"
    a.cell.cyt2$rep <- 2
    a.cell.nuc2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "nucleus"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc2$loc <- "nucleus"
    a.cell.nuc2$rep <- 2
    
    comb <- rbind(a.cell.cyt1,a.cell.cyt2,a.cell.nuc1,a.cell.nuc2)
    transExpr <- as.data.frame(group_by(comb,gene_id,transcript_id) %.% summarise(readsPerLen=sum(readsPerLen)))
    colnames(transExpr) <- c("transcript_id", "gene_id", "readsPerLen")
    transExpr$readsPerLen_tieBreaker <- transExpr$readsPerLen + runif(seq_along(transExpr$readsPerLen))/(10^9)
    gene.df <- as.data.frame(group_by(transExpr, gene_id) %.% filter(readsPerLen_tieBreaker == max(readsPerLen_tieBreaker))) 
    cellTranscripts <- gene.df$transcript_id
    cTrans <- comb[which(comb$transcript_id %in% cellTranscripts),]
    cTrans$cell <- cell                       
    
    df.together <- rbind(df.together,cTrans)
    
  }
  exportAsTable(file=getFullPath("data/rpkmFromBam-TopTransCellType-Reads.tab"),df=df.together)
}