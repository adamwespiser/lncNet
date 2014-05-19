
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
                         "rpkmFromBamFile", "bare", "readLength","rfbTrans")
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
  1
}

processCellsMaxTransExpr_ByExon <- function(){
  annot.df <- getRpkmFromBamDataForOneCellByExon()
  annot.df <- annot.df[which(annot.df$rnaExtract == "longPolyA"),]
  annot.df <- annot.df[-which(annot.df$cell == "H1-hESC"),]
  trans <- annot.df$rfbTrans
  exons <- annot.df$rpkmFromBamFile
  annot.df$rep <- ifelse(annot.df$replicate >2,annot.df$replicate -2,annot.df$replicate )
#  annot.df$transFullReads<- gsub(x=annot.df$rfbGene,pattern="genes",replacement="transFullReads")
 # transFullReads <- annot.df$transFullReads 
  sapply(seq_along(exons), function(x)convertExonsToTrans(exonFile=exons[x],transFile=trans[x]))
  df.together <- data.frame()
  for ( cell in unique(annot.df$cell)){
    print(cell)
    a.cell <- annot.df[which(annot.df$cell == cell),]
    a.cell.cyt1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "cytosol"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt1$localization <- "cytosol"
    a.cell.cyt1$replicate <- 1
    a.cell.nuc1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "nucleus"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc1$localization <- "nucleus"
    a.cell.nuc1$replicate <- 1
    a.cell.cyt2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "cytosol"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt2$localization <- "cytosol"
    a.cell.cyt2$replicate <- 2
    a.cell.nuc2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "nucleus"),"rfbTrans"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc2$localization <- "nucleus"
    a.cell.nuc2$replicate <- 2
    
    comb <- rbind(a.cell.cyt1,a.cell.cyt2,a.cell.nuc1,a.cell.nuc2)
    transExpr <- as.data.frame(group_by(comb,gene_id,transcript_id) %.% summarise(readsPerLen=sum(readsPerLen)))
    #colnames(transExpr) <- c("transcript_id", "gene_id", "readsPerLen")
    transExpr$readsPerLen_tieBreaker <- transExpr$readsPerLen + runif(seq_along(transExpr$readsPerLen))/(10^9)
    gene.df <- as.data.frame(group_by(transExpr, gene_id) %.% filter(readsPerLen_tieBreaker == max(readsPerLen_tieBreaker))) 
    cellTranscripts <- gene.df$transcript_id
    cTrans <- comb[which(comb$transcript_id %in% cellTranscripts),]
    cTrans$cell <- cell                       
    cTrans$reads <- as.numeric(cTrans$reads)
    cTrans$length <- as.numeric(cTrans$length)
    cTransRPKM <- as.data.frame(group_by(cTrans,localization,replicate) %.% mutate(millReads = sum(reads,na.rm=TRUE)/(10^6),
                                                                                   RPKM = (reads * 10^9)/(length * sum(reads,na.rm=TRUE)) ))
    
    df.together <- rbind(df.together,cTransRPKM)
    
  }
  exportAsTable(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM.tab"),df=df.together)
}

getDataTotalReadsBtwnReps_rpkmFromBamTopTrans_ByExon <- function(){
  df.together <- read.csv(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM.tab"),sep="\t")
  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  df.together$region <- "other"
  df.together[which(df.together$gene_id %in% lnc),"region"] <- "lnc"
  df.together[which(df.together$gene_id %in% pc),"region"] <- "mRNA"
  df.together$gene_type <- df.together$region
  df.together$rnaExtract = "longPolyA"

  df.together$RPKM <- as.numeric(df.together$RPKM)
  df.together$isSpikeIn <- 0
  df.together[grep(pattern="ERCC",df.together$gene_id),"isSpikeIn"] <- 1
  
  report.df  <- as.data.frame(group_by(df.together,cell,localization,replicate) %.%
                                summarise(length(gene_id),
                                          mean(RPKM),
                                          sum(RPKM),
                                          sum(RPKM > 0)))
  report.df$experiment <- paste(ifelse(report.df$localization == "cytosol", "cyt", "nuc"),report.df$replicate,sep=".")
  colnames(report.df) <- c("cell", "localization", "replicate", "genesFound", "meanRPKM", 
                           "sumRPKM","genesExpressed", "experiment")
  exportAsTable(df=report.df, file = getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-REPORT.tab"))
  
  df.together <- as.data.frame(group_by(df.together, cell, localization,rnaExtract,replicate) %.% 
                                 mutate(RPKM_80norm = apply80norm(RPKM) * 1000000))
  
  #  group_by(df.together, cell, localization,rnaExtract,replicate) %.% summarise(mean(RPKM_80norm/transTotalRPKM, na.rm=TRUE))
  
  #exportAsTable(file=getFullPath("/data/rpkmFromBamTopTransAllCells.tab"), df=df.together)
  df.together$gene_type <- df.together$region
  df.abbrev <- df.together[ c("region","replicate", "gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn", "RPKM_80norm","RPKM","reads")]
  
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
  
  
  #df.cytNuc <- rbind(df.cyt,df.nuc)
  #df.cytNuc[which(df.cytNuc$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc[which(df.cytNuc$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc <- merge(df.cyt,df.nuc,by=c("gene_id","cell","variable"),suffixes=c(".cyt",".nuc"))
  
  
  exportAsTable(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-proc.tab"), df=df.cytNuc)
  
}


plotRatiosTopTrans <- function(){
  
  df.cytNuc <- read.csv(sep="\t",file=getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-proc.tab"))
  df.cytNuc$cytFracPseudo <- with(df.cytNuc, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  df.cytNuc$cytFrac <- with(df.cytNuc, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  df.cytNuc.pos <- df.cytNuc[which(df.cytNuc$value.ave.cyt != 0 & df.cytNuc$value.ave.nuc != 0),]
  
  
  df.lpa.ratio.rpkm <- df.cytNuc.pos[which(df.cytNuc.pos$variable =="RPKM"),]
  df.lpa.ratio.rpkm80 <- df.cytNuc.pos[which(df.cytNuc.pos$variable =="RPKM_80norm"),]
  
  if(!file.exists(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/"))){
    dir.create(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/"))
  }
  
  
  ggplot(df.lpa.ratio.rpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkmPseudo-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cyt/(nuc + cyt)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkm-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm, aes(x=cytFracPseudo,fill=factor(region.cyt)))+
    geom_bar(position="dodge") + theme_bw() + thisTheme + 
    facet_grid(cell~.)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkmPseudo-bars-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm, aes(x=cytFrac,,fill=factor(region.cyt)))+
    geom_bar(position="dodge") + theme_bw() + thisTheme + 
    facet_grid(cell~.)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cyt/(nuc + cyt)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkm-bars-cells.png"), height=12,width=5)
  
  
  #RPKM80 norm
  ggplot(df.lpa.ratio.rpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM80: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkm80Pseudo-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM80: cyt/(nuc + cyt)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/rpkm80-cells.png"), height=12,width=5)
  
  
}




