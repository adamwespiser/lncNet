



plotMultiVsUniqReads <- function(){
  
  
  rfbMulti <- read.csv(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-proc.tab"), stringsAsFactors=FALSE, sep ="\t")
  #uniqOutfile=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM-proc.tab")
  rfbUniq <- read.csv(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM-proc.tab"), stringsAsFactors=FALSE, sep ="\t")
  
  comb <- merge(x=rfbMulti, y=rfbUniq, by=c("gene_id","cell","variable"),suffixes=c(".multi",".uniq"))
  
  
  ggplot(comb[which(comb$variable == "RPKM"),], aes(x=log10(2*value.ave.nuc.multi + 2*value.ave.cyt.multi),
                                                    y=log10(2*value.ave.nuc.uniq + 2*value.ave.cyt.uniq),
                                                    color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM -> uniq vs. multi mapped reads\n")+
    xlab("multi reads over cyt+nuc rep1+2")+
    ylab("unique reads over cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompare/uniqVsMulti-combined.png"),height=8,width=8)
  

  ggplot(comb[which(comb$variable == "TPM"),], aes(x=log10(2*value.ave.nuc.multi + 2*value.ave.cyt.multi),
                                                     y=log10(2*value.ave.nuc.uniq + 2*value.ave.cyt.uniq),
                                                     color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM -> uniq vs. multi mapped reads\n")+
    xlab("multi reads over cyt+nuc rep1+2")+
    ylab("unique reads over cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() +
    facet_wrap(~cell)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompare/tpm-uniqVsMulti-byCell.png"),height=8,width=8)
  
  
  ggplot(comb[which(comb$variable == "readsPerKb"),], aes(x=log10(2*value.ave.nuc.multi + 2*value.ave.cyt.multi),
                                                     y=log10(2*value.ave.nuc.uniq + 2*value.ave.cyt.uniq),
                                                     color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM -> uniq vs. multi mapped reads\n")+
    xlab("multi reads over cyt+nuc rep1+2")+
    ylab("unique reads over cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() +
    facet_wrap(~cell)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompare/readsPerKb-uniqVsMulti-byCell.png"),height=8,width=8)
  ggplot(comb[which(comb$variable == "RPKM"),], aes(x=log10(2*value.ave.nuc.multi + 2*value.ave.cyt.multi),
                                                    y=log10(2*value.ave.nuc.uniq + 2*value.ave.cyt.uniq),
                                                    color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM -> uniq vs. multi mapped reads\n")+
    xlab("multi reads over cyt+nuc rep1+2")+
    ylab("unique reads over cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompare/uniqVsMulti-combined.png"),height=8,width=8)
  
  
  #comb.reads <- comb[which(comb$variable == "reads"),]
  #comb.reads$multi <- 2*comb.reads$value.ave.nuc.multi + 2*comb.reads$value.ave.cyt.multi
  #comb.reads$uniq <- 2*comb.reads$value.ave.nuc.uniq + 2*comb.reads$value.ave.cyt.uniq
  #annom <- comb.reads[which(comb.reads$uniq > comb.reads$multi),]
  
  #subset(rfbMulti, gene_id == "ENSG00000001626.10" & cell == "MCF-7")[c("variable","value.ave.nuc","value.ave.cyt")]
  #subset(rfbUniq, gene_id == "ENSG00000001626.10" & cell == "MCF-7")[c("variable","value.ave.nuc","value.ave.cyt")]
 
}

rpkmFromBamVsRSEM <- function(){
  rfbUniq <- read.csv(file=getFullPath("data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM-proc.tab"), stringsAsFactors=FALSE, sep ="\t")
  rfbUniq$cytFracPseudo <- with(rfbUniq, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  rfbUniq$cytFrac <- with(rfbUniq, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  rfbUniq$rnaExtract <- rfbUniq$rnaExtract.cyt
 
  RSEM <-  read.csv(file=getFullPath("/data/rsemCapData-v2-lpa-proc.tab.merge"), stringsAsFactors=FALSE, sep ="\t")
  RSEM$variable <- gsub(x=RSEM$variable, pattern="FPKM",replacement="RPKM")
  RSEM$cytFracPseudo <- with(RSEM, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  RSEM$cytFrac <- with(RSEM, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  comb <- merge(x=RSEM, y=rfbUniq, by=c("gene_id","cell","rnaExtract.nuc","variable"),suffixes=c(".rsem",".rfb"))
  
  
  ggplot(comb[which(comb$variable == "RPKM"),], aes(x=log10(2*value.ave.nuc.rsem + 2*value.ave.cyt.rsem),
                                                    y=log10(2*value.ave.nuc.rfb + 2*value.ave.cyt.rfb),
                                                    color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM uniq vs. RSEM multi mapped reads\n")+
    xlab("RSEM: cyt+nuc rep1+2")+
    ylab("RPKMfromBAM cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompareRSEM/combined.png"),height=8,width=8)
  
  
  ggplot(comb[which(comb$variable == "TPM"),], aes(x=log10(2*value.ave.nuc.rsem + 2*value.ave.cyt.rsem),
                                                   y=log10(2*value.ave.nuc.rfb + 2*value.ave.cyt.rfb),
                                                   color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM vs. RSEM\n")+
    xlab("RSEM: cyt+nuc rep1+2")+
    ylab("RPKMfromBAM cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() +
    facet_wrap(~cell)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompareRSEM/tpm-byCell.png"),height=8,width=8)
  
  
  ggplot(comb[which(comb$variable == "RPKM"),], aes(x=cytFrac.rsem,
                                                   y=cytFrac.rfb,
                                                   color=cell)) +
    geom_density2d() + geom_point(size=1)+
    ggtitle("RPKMfromBAM vs. RSEM\n")+
    xlab("RSEM: cytosolic Frac.(RPKM)")+
    ylab("RPKMfromBAM cytosolic Frac.(RPKM)")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() +
    facet_grid(cell~.)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompareRSEM/rpkm-cytFrac-byCell.png"),height=12,width=8)
  
  
  # to update these plots, re-run RSEM-datainput scripts, all files should update...
  
  
  
  
  starReads <- read.csv(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/STAR-mappedStats"),sep="\t",stringsAsFactors=FALSE)
  starReads$rep <- ifelse(starReads$replicate > 2, starReads$replicate -2 , starReads$replicate)
  starReads.lpa <- starReads[which(starReads$rnaExtract == "longPolyA"),] 
  starReads.lpa.nuc <- starReads[which(starReads$rnaExtract == "longPolyA" & starReads$localization == "nucleus"),] 
  starReads.lpa.cyt<- starReads[which(starReads$rnaExtract == "longPolyA" & starReads$localization == "cytosol"),] 
  rfbMappedReads <- as.data.frame(group_by(starReads,cell,localization, rnaExtract) %.% summarise(mappedReads = sum(uniqMapped)))
  rfbMappedReads.lpa <- rfbMappedReads[which(rfbMappedReads$rnaExtract == "longPolyA"),]
  rfbMappedReads.lpa$rnaExtract <- NULL
  melt(rfbMappedReads, id.var="cell")
  
  rfbUniqCorr <- melt(group_by(rfbUniq[which(rfbUniq$variable == "RPKM"),],cell) %.% 
                        summarise(nucleus =cor(value.rep1.nuc,value.rep2.nuc,method="spearman"),
                                  cytosol =cor(value.rep1.cyt,value.rep2.cyt,method="spearman")),id.var="cell")  
  colnames(rfbUniqCorr) <- c("cell","localization","spearmanCorr")
  
  rfbUniqMap <- merge(rfbMappedReads.lpa,rfbUniqCorr, by=c("cell","localization"))
  
  ggplot(rfbUniqMap, aes(x=mappedReads,y=spearmanCorr,label=paste(cell,localization,sep="  "),
                         color=paste(cell,localization),size=6)) + 
    geom_point()+geom_text()+
    theme_bw() + 
    xlab("Number of STAR mapped unique Reads") + 
    ylab("Spearman Correlation between replicate")+
    ggtitle("STAR + rpkmFromBam : replicate correlation vs. mappedReads")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompareRSEM/rpkmFromBam-corrVsReadCount.png"),height=8,width=8)
  
  
  
  
  rsemDistro.file <- getFullPath("/data/rsemCapData-v2-readDistro.tab")
  rsemReads <- read.csv(file=rsemDistro.file, stringsAsFactors=FALSE,sep="\t")
  rsemReads <- unique(rsemReads[c("cell","localization","replicate","rnaExtract","mapped")])
  rsemMappedReads <- as.data.frame(group_by(rsemReads,cell,localization, rnaExtract) %.% summarise(mappedReads = sum(mapped)))
  
  
  rsemMappedReads.lpa <- rsemMappedReads[which(rsemMappedReads$rnaExtract == "longPolyA"),]
  rsemMappedReads.lpa$rnaExtract <- NULL
  
  
  RSEMCorr <- melt(group_by(RSEM[which(RSEM$variable == "RPKM"),],cell) %.% 
                        summarise(nucleus =cor(value.rep1.nuc,value.rep2.nuc,method="spearman"),
                                  cytosol =cor(value.rep1.cyt,value.rep2.cyt,method="spearman")),id.var="cell")  
  colnames(RSEMCorr) <- c("cell","localization","spearmanCorr")
  RSEMplot <- merge(RSEMCorr, rsemMappedReads.lpa, by = c("cell","localization"))
  ggplot(RSEMplot, aes(x=mappedReads,y=spearmanCorr,label=paste(cell,localization,sep="  "),
                         color=paste(cell,localization),size=6)) + 
    geom_point()+geom_text()+
    theme_bw() + 
    xlab("Number of RSEM mapped Reads") + 
    ylab("Spearman Correlation between replicate")+
    ggtitle("RSEM(bwa mapping) : replicate correlation vs. mappedReads")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompareRSEM/rsem-corrVsReadCount.png"),height=8,width=8)
  
  
  RSEMCorr.pearson <- melt(group_by(RSEM[which(RSEM$variable == "RPKM"),],cell) %.% 
                     summarise(nucleus =cor(value.rep1.nuc,value.rep2.nuc,method="pearson"),
                               cytosol =cor(value.rep1.cyt,value.rep2.cyt,method="pearson")),id.var="cell")  
  
}



