





plotDifferenceBetweenRepsRSEM <- function(file=getFullPath("/data/rsemCapData-lpa-proc.tab") )  {
  
  df.cytNuc <- read.csv(file=file, stringsAsFactors=FALSE, sep ="\t")
  #df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.fpkm <- df.cytNuc[which(df.cytNuc$variable =="FPKM"),]
  #df.cytNuc.fpkm[which(df.cytNuc.fpkm$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.fpkm[which(df.cytNuc.fpkm$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.fpkm80<- df.cytNuc[which(df.cytNuc$variable =="FPKM_80norm"),]
  #transcriptTotalConc
  df.cytNuc.tpm<- df.cytNuc[which(df.cytNuc$variable =="TPM"),]
  
  
  
  
  # rep vs. rep LOG
  ggplot(df.cytNuc.tpm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nTPM \nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("TPM") + ylab("TPM")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/tpm-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.fpkm80, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nFPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM 80 norm.") + ylab("FPKM 80 norm")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/fpkm80-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nFPKM\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM") + ylab("FPKM ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/fpkm-vsReps.png"), height=12,width=5)
  
  
  # rep vs. rep LOG
  ggplot(df.cytNuc.tpm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nTPM n\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("TPM") + ylab("TPM")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/tpm-log-vsReps.png"), height=12,width=5)
  
    
  ggplot(df.cytNuc.fpkm80, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nFPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM 80 norm.") + ylab("FPKM 80 norm")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/fpkm80-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RSEM\nFPKM \nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM") + ylab("FPKM ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/fpkm-log-vsReps.png"), height=12,width=5)
  
  
  
  #mRNA  only 
#   ggplot(df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$region == "mRNA"),], 
#          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + 
#     geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RSEM\nFPKM normalized by Spike In\nLongPolyA only\nmRNA only")+ 
#     geom_abline(slope=1,intercept=0)+ 
#     xlab("RPKM spike norm.") + ylab("RPKM spike norm")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkmSpikeIn-mRNAOnly-log-vsReps.png"), height=12,width=5)
#   
#   ggplot(df.cytNuc.fpkm[which(df.cytNuc.fpkm$region == "mRNA"),],
#          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RSEM\nRPKM of Spike Ins\nmRNA only\nLongPolyA only")+ 
#     geom_abline(slope=1,intercept=0)+ xlab("RPKM") + ylab("RPKM ")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkm-mRNAOnly-log-vsReps.png"), height=12,width=5)
#   
#   #lncRNA
#   ggplot(df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$region == "lncRNA"),],
#          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RSEM\nRPKM normalized by Spike In\nlncRNA only\nLongPolyA only")+ 
#     geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm.") + ylab("RPKM spike norm")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkmSpikeIn-lncOnly-log-vsReps.png"), height=12,width=5)
#   
#   ggplot(df.cytNuc.fpkm[which(df.cytNuc.fpkm$region == "lncRNA"),], 
#          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RSEM\nRPKM of Spike Ins\nlncRNA only\nLongPolyA only")+ 
#     geom_abline(slope=1,intercept=0)+ xlab("RPKM") + ylab("RPKM ")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkm-lncOnly-log-vsReps.png"), height=12,width=5)
#   
  
  
  m.df <- as.data.frame(group_by(df.cytNuc, cell, localization,variable,region) %.%
                          
                          summarize(sum.rep1 = sum(value.rep1),
                                    sum.rep2 = sum(value.rep2),
                                    expr.rep1 = sum(value.rep1 > 0 ),
                                    expr.rep2 = sum(value.rep2 > 0 )))
  
  colnames(m.df) <- c("cell", "localization", "measure", "region", "sum.rep1", "sum.rep2", "expr.rep1", "expr.rep2")
  melt.df <- melt(m.df, id.vars=c("cell", "localization", "measure","region"))
  m.df$frac.rep1 = with(m.df, (sum.rep1)/(sum.rep1 + sum.rep2))
  
  
  ggplot(melt.df[which(melt.df$variable %in% c("expr.rep1", "expr.rep2") & melt.df$measure == "TPM"),], aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + 
    facet_grid(cell ~localization) + 
    ggtitle("RSEM\nnumber of expressed genes")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/expr-cytNuc.png"), height=12,width=5)
  # localization vs. cell 
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) +  
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(localization~cell) + 
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme + ggtitle("RSEM\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/readCount-cytNuc-region.png"), height=5,width=12)
  
  # combined
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle("RSEM\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/readCount-cytNuc-all-combined.png"), height=5,width=10)
  
  ggplot(m.df, aes(x=measure,y=frac.rep1)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle("RSEM\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/readCount-cytNuc-all-combined-join.png"), height=5,width=10)
  
  
  #seperate by localization only
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(~localization) +
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle("RSEM\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/readCount-cytNuc-combined.png"), height=6,width=10)
  
  
  
  
}
