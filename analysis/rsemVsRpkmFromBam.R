







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
  
  
  ggplot(comb[which(comb$variable == "reads"),], aes(x=log10(2*value.ave.nuc.multi + 2*value.ave.cyt.multi),
                                                     y=log10(2*value.ave.nuc.uniq + 2*value.ave.cyt.uniq),
                                                     color=cell)) +
    geom_point()+
    ggtitle("RPKMfromBAM -> uniq vs. multi mapped reads\n")+
    xlab("multi reads over cyt+nuc rep1+2")+
    ylab("unique reads over cyt+nuc rep1+2")+
    geom_abline(slope=1,intercept=0) +
    theme_bw() +
    facet_wrap(~cell)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBamCompare/reads-uniqVsMulti-byCell.png"),height=8,width=8)
  
  
  #comb.reads <- comb[which(comb$variable == "reads"),]
  #comb.reads$multi <- 2*comb.reads$value.ave.nuc.multi + 2*comb.reads$value.ave.cyt.multi
  #comb.reads$uniq <- 2*comb.reads$value.ave.nuc.uniq + 2*comb.reads$value.ave.cyt.uniq
  #annom <- comb.reads[which(comb.reads$uniq > comb.reads$multi),]
  
  #subset(rfbMulti, gene_id == "ENSG00000001626.10" & cell == "MCF-7")[c("variable","value.ave.nuc","value.ave.cyt")]
  #subset(rfbUniq, gene_id == "ENSG00000001626.10" & cell == "MCF-7")[c("variable","value.ave.nuc","value.ave.cyt")]
  
  
}