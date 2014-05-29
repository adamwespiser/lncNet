


doit <- function(){
  processCellsMaxTransExpr_ByExon()
  getDataTotalReadsBtwnReps_rpkmFromBamTopTrans_ByExon()
}

doitUniqReads <- function(){
  getRpkmFromBamDataForOneCellByExon(suffix=".uniq.star_sort.transByExon.gtf",transSuffix=".uniq.star_sort.trans.gtf",writeCopyScript=TRUE)
  system("~/sandbox/rpkmFromBamExonFetch")
  
  uniqOutfile=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM.tab")
  uniqReportFile = getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM-REPORT.tab")
  processCellsMaxTransExpr_ByExon(suffix=".uniq.star_sort.transByExon.gtf",transSuffix=".uniq.star_sort.trans.gtf",outfile=uniqOutfile)
  procFile = getDataTotalReadsBtwnReps_rpkmFromBamTopTrans_ByExon(infile=uniqOutfile,reportFile=uniqReportFile)
  plotRatiosTopTrans(infile=procFile,
                     outdir = getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon_uniqRead/"),
                     plotMsg="unique reads mapped by STAR")
  
}

getRpkmFromBamDataForOneCellByExon <- function( filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",writeCopyScript=FALSE,suffix=".transByExon.gtf",transSuffix){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df$readLength <- as.numeric(gsub(sapply(strsplit(df$readType, "x"), function(x)x[2]),pattern="D",replacement=""))
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  #/project/umw_zhiping_weng/wespisea/rna-seq//starSpikeIn/
  df.comb$remote <- file.path(rnaseqdir,"starSpikeIn/",paste0(df.comb$bare,suffix))
  
  if(writeCopyScript){
    o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, " /home/wespisea/data/rpkmFromBam/",paste0(df.comb$bare,suffix))
    write(o1,file="~/sandbox/rpkmFromBamExonFetch")
  }
  df.comb$rpkmFromBamFile <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,suffix)
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"rpkmFromBamFile", "bare","read1.readLength")]
  df.comb$rfbTrans <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,transSuffix)
  colnames(df.comb) <- c("localization", "cell", "rnaExtract", "replicate", 
                         "rpkmFromBamFile", "bare", "readLength","rfbTrans")
  df.comb
}


convertExonsToTrans <- function(transFile,exonFile){
  if(!file.exists(exonFile)){
    return(0)
  }
  
  print(paste(transFile,exonFile,sep="\t"))
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
  trans.df$readsPerKb <- with(trans.df, reads/(length/1000))
  
  spike.df <- getSpikeInDf()
  
  
  df.spike <- trans.df[grep(pattern="ERCC",trans.df$gene_id),]
  if(dim(df.spike)[1] == 0){
    trans.df$concBySpikeIn <- NA
    
  } else {
    df.spike$readsPerKb <- df.spike$reads / (df.spike$length/1000)
    #df.spike$RPKM <- df.spike$readsPerKb / millionsOfReads
    spike <- merge(df.spike,spike.df, by="gene_id")
    #s.lm <- glm(Pool14nmol.ul ~ reads + length, data=spike,family="poisson")
    s.lm <- glm(Pool14nmol.ul ~ readsPerKb + 0, data=spike)
    trans.df$concBySpikeIn <- predict(s.lm, newdata=trans.df)
    
  }
  
  
  
  exportAsTable(df=trans.df  ,file=transFile)
  1
}

processCellsMaxTransExpr_ByExon <- function(suffix=".transByExon.gtf",
                                            transSuffix=".transFromExon.gtf",
                                            skipCells=c("none"),
                                            outfile=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM.tab")){
  annot.df <- getRpkmFromBamDataForOneCellByExon(suffix=suffix,transSuffix=transSuffix)
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
    readLen.vec <-  annot.df[which(annot.df$cell == cell),"readLength"]
    if (FALSE == (sum(readLen.vec[1] == readLen.vec) == length(readLen.vec))){
      print(paste("READ LENGTHS DIFFER :-( ", cell, "problematic"))
    }
    a.cell <- annot.df[which(annot.df$cell == cell),]
    cyt1.file <- a.cell[which(a.cell$rep == 1 & a.cell$localization == "cytosol"),"rfbTrans"]
    cyt2.file <- a.cell[which(a.cell$rep == 2 & a.cell$localization == "cytosol"),"rfbTrans"]
    nuc1.file <- a.cell[which(a.cell$rep == 1 & a.cell$localization == "nucleus"),"rfbTrans"]
    nuc2.file <- a.cell[which(a.cell$rep == 2 & a.cell$localization == "nucleus"),"rfbTrans"]
    
    if(file.exists(cyt1.file) &&  file.exists(cyt2.file) && file.exists(nuc1.file) && file.exists(nuc2.file)){
      
      #cyt rep 1
      a.cell.cyt1 <- read.csv(file=cyt1.file,sep="\t",stringsAsFactors=FALSE)
      a.cell.cyt1$localization <- "cytosol"
      a.cell.cyt1$replicate <- 1
      # cyt rep 2
      a.cell.nuc1 <- read.csv(file=cyt2.file,sep="\t",stringsAsFactors=FALSE)
      a.cell.nuc1$localization <- "cytosol"
      a.cell.nuc1$replicate <- 2
      # nuc rep1 
      a.cell.cyt2 <- read.csv(file=nuc1.file,sep="\t",stringsAsFactors=FALSE)
      a.cell.cyt2$localization <- "nucleus"
      a.cell.cyt2$replicate <- 1
      #nuc rep 2
      a.cell.nuc2 <- read.csv(file=nuc2.file,sep="\t",stringsAsFactors=FALSE)
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
      cTrans$readLength = readLen.vec[1]
      cTransRPKM <- as.data.frame(group_by(cTrans,localization,replicate) %.% mutate(millReads = sum(reads,na.rm=TRUE)/(10^6),
                                                                                     RPKM = (reads * 10^9)/(length * sum(reads,na.rm=TRUE)),
                                                                                     TPM  = (reads * readLength * 10^6)/(length * (sum(reads * readLength/length))) ))
      
      df.together <- rbind(df.together,cTransRPKM)
    }
  }
  exportAsTable(file=outfile,df=df.together)
}

getDataTotalReadsBtwnReps_rpkmFromBamTopTrans_ByExon <- function(infile=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM.tab"),
                                                                 reportFile =getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-REPORT.tab") ){
  df.together <- read.csv(file=infile,sep="\t")
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
                                          sum(RPKM > 0),
                                          sum(isSpikeIn)))
  
  report.df$experiment <- paste(ifelse(report.df$localization == "cytosol", "cyt", "nuc"),report.df$replicate,sep=".")
  colnames(report.df) <- c("cell", "localization", "replicate", "genesFound", "meanRPKM", 
                           "sumRPKM","genesExpressed","spikeInDetected" ,"experiment")
  exportAsTable(df=report.df, file = reportFile)
  
  
  reportTPM.df  <- as.data.frame(group_by(df.together,cell,localization,replicate) %.%
                                   summarise(length(gene_id),
                                             mean(RPKM),
                                             sum(RPKM),
                                             sum(RPKM > 0),
                                             sum(isSpikeIn)))
  
  reportTPM.df$experiment <- paste(ifelse(reportTPM.df$localization == "cytosol", "cyt", "nuc"),reportTPM.df$replicate,sep=".")
  colnames(reportTPM.df) <- c("cell", "localization", "replicate", "genesFound", "meanTPM", 
                              "sumTPM","genesExpressed","spikeInDetected" ,"experiment")
  exportAsTable(df=report.df, file = getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-TPM-REPORT.tab"))
  
  
  
  df.together <- as.data.frame(group_by(df.together, cell, localization,rnaExtract,replicate) %.% 
                                 mutate(RPKM_80norm = apply80norm(RPKM) * 1000000))
  
  #  group_by(df.together, cell, localization,rnaExtract,replicate) %.% summarise(mean(RPKM_80norm/transTotalRPKM, na.rm=TRUE))
  
  #exportAsTable(file=getFullPath("/data/rpkmFromBamTopTransAllCells.tab"), df=df.together)
  df.together$gene_type <- df.together$region
  df.abbrev <- df.together[ c("region","replicate", "gene_id","gene_type", "localization","rnaExtract","cell", "isSpikeIn", "RPKM_80norm","RPKM","reads","concBySpikeIn","TPM")]
  
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
  df.cytNuc.rbind <- rbind(df.cyt,df.nuc)
  
  in.vec <- unlist(strsplit(infile,"\\."))
  outfile = paste0(in.vec[1],"-proc.",in.vec[2])
  
  
  exportAsTable(file=outfile, df=df.cytNuc)
  exportAsTable(file=paste0(outfile,".rbind"), df=df.cytNuc.rbind)
  
  outfile
}


plotRatiosTopTrans <- function(infile,outdir = getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamTopTrans/cytFracByExon/"),
                                                           plotMsg=" "){
  
  df.cytNuc <- read.csv(sep="\t",file=infile)
  df.cytNuc$cytFracPseudo <- with(df.cytNuc, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  df.cytNuc$cytFrac <- with(df.cytNuc, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  df.cytNuc.pos <- df.cytNuc[which(df.cytNuc$value.ave.cyt != 0 & df.cytNuc$value.ave.nuc != 0),]
  
  
  df.lpa.ratio.rpkm <- df.cytNuc.pos[which(df.cytNuc.pos$variable =="RPKM"),]
  df.lpa.ratio.tpm <- df.cytNuc.pos[which(df.cytNuc.pos$variable =="TPM"),]
  
  df.lpa.ratio.rpkm80 <- df.cytNuc.pos[which(df.cytNuc.pos$variable =="RPKM_80norm"),]
  
  if(!file.exists(paste0(outdir,"/"))){
    dir.create(outdir)
  }
  
  
  ggplot(df.lpa.ratio.rpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cytPseudo/(nucPseudo + cytPseudo)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkmPseudo-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste0("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cyt/(nuc + cyt)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkm-cells.png"), height=12,width=5)
  
  
  #TPM 
  ggplot(df.lpa.ratio.tpm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nTPM: cytPseudo/(nucPseudo + cytPseudo)\n",plotMsg))
  ggsave(paste0(outdir,"/tpmPseudo-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.tpm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nTPM: cyt/(nuc + cyt)\n",plotMsg))
  ggsave(paste0(outdir,"/tpm-cells.png"), height=12,width=5)
  
  
  
  ggplot(df.lpa.ratio.rpkm, aes(x=cytFracPseudo,fill=factor(region.cyt)))+
    geom_bar(position="dodge") + theme_bw() + thisTheme + 
    facet_grid(cell~.)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cytPseudo/(nucPseudo + cytPseudo)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkmPseudo-bars-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm, aes(x=cytFrac,,fill=factor(region.cyt)))+
    geom_bar(position="dodge") + theme_bw() + thisTheme + 
    facet_grid(cell~.)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM: cyt/(nuc + cyt)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkm-bars-cells.png"), height=12,width=5)
  
  
  #RPKM80 norm
  ggplot(df.lpa.ratio.rpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM80: cytPseudo/(nucPseudo + cytPseudo)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkm80Pseudo-cells.png"), height=12,width=5)
  
  ggplot(df.lpa.ratio.rpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle(paste("RPKMfromBAM Top Trans Count By Exon\nFraction of Cytosolic RNA-seq expr\nRPKM80: cyt/(nuc + cyt)\n",plotMsg))
  ggsave(paste0(outdir,"/rpkm80-cells.png"), height=12,width=5)
  
  
  df.cytNuc.rbind <- read.csv(sep="\t",file=paste0(infile,".rbind"))
  
  m.df <- as.data.frame(group_by(df.cytNuc.rbind , cell, localization,variable,region) %.%
                          
                          summarize(sum.rep1 = sum(value.rep1),
                                    sum.rep2 = sum(value.rep2),
                                    expr.rep1 = sum(value.rep1 > 0 ),
                                    expr.rep2 = sum(value.rep2 > 0 )))
  
  colnames(m.df) <- c("cell", "localization", "measure", "region", "sum.rep1", "sum.rep2", "expr.rep1", "expr.rep2")
  melt.df <- melt(m.df, id.vars=c("cell", "localization", "measure","region"))
  m.df$frac.rep1 = with(m.df, (sum.rep1)/(sum.rep1 + sum.rep2))
  
  
  ggplot(melt.df[which(melt.df$variable %in% c("expr.rep1", "expr.rep2") & melt.df$measure == "RPKM"),], aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + 
    facet_grid(cell ~localization) + xlab("replicates") + ylab("count") +
    ggtitle(paste("RPKMfromBam Top Trans Count by Exon\nnumber of expressed genes\n",plotMsg))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"/expr-cytNuc.png"), height=12,width=5)
  # localization vs. cell 
  # "RPKM_80norm","RPKM","reads","concBySpikeIn","spikeIn_norm
  
  
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) +  
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(localization~cell) + 
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","concBySpikeIn","TPM"),
                     labels=c("reads","RPKM"  ,"RPKM_80","conc.","TPM"))+ ylim(0,1) +
    thisTheme + 
    ggtitle(paste("RPKMfromBam Top Trans Count by Exon\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/conc/TPM\nfrac.rep1=(rep1)/(rep1 + rep2)\n",plotMsg))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"/readCount-cytNuc-region.png"), height=5,width=12)
  
  # combined
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","concBySpikeIn","TPM"),
                     labels=c("reads","RPKM"  ,"RPKM_80","conc.","TPM"))+ ylim(0,1) +
    thisTheme2 + 
    ggtitle(paste("RPKMfromBam Top Trans Count by Exon\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMspikeIn/Conc/TPM\n frac.rep1=(rep1)/(rep1 + rep2)\n",plotMsg))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"/readCount-cytNuc-all-combined.png"), height=5,width=10)
  
  ggplot(m.df, aes(x=measure,y=frac.rep1)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","concBySpikeIn","TPM"),
                     labels=c("reads","RPKM"  ,"RPKM_80","conc.","TPM"))+ ylim(0,1) +
    thisTheme2 + 
    ggtitle(paste("RPKMfromBam Top Trans Count by Exon\nfraction of cytosol & nucleus\nreads/RPKM/RPKM80/RPKMspikeIn/Conc/TPM\nfrac.rep1=(rep1)/(rep1 + rep2)\n",plotMsg))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"/readCount-cytNuc-all-combined-join.png"), height=5,width=10)
  
  
  #seperate by localization only
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(~localization) +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","concBySpikeIn","TPM"),
                     labels=c("reads","RPKM"  ,"RPKM_80","conc.","TPM"))+ ylim(0,1) +
    thisTheme2 + ggtitle(paste("RPKMfromBam Top Trans Count by Exon\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMspikeIn/Conc/TPM\nfrac.rep1=(rep1)/(rep1 + rep2)\n",plotMsg))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"/readCount-cytNuc-combined.png"), height=6,width=10)
  
}




