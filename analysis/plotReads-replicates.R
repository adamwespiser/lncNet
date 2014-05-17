


plotSpikeInOnly <- function(file=getFullPath("./data/spikeIn-cytNucPapPam.tab")){
  df.cytNuc <- read.csv(file=file, sep="\t", stringsAsFactors=FALSE)
  # determine the relationship between fraction of reads(expected value = 1/2) and nucleotides
  df.cytNuc.reads <- df.cytNuc[which(df.cytNuc$variable == "reads"),]
  
  
  
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(rep1.frac.pseudo ~ nt, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(nt, rep1.frac.pseudo)) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme +
    geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("length of spikeIn- in nt") +
    ggtitle("long Poly A pulldown \n Rep 1 Fraction of Total Reads\nFrac = (rep1 + pCount)/(rep1 + rep2 + 2*pCount)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/fracError-by-nt.pdf"), height=12,width=5)
  # determine the relationship between fraction of reads(expected value = 1/2) and gc concent
  #df.nuc.reads <- df.nuc[which(df.nuc$variable == "reads"),]
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(rep1.frac.pseudo ~ percentGC, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(x=percentGC, rep1.frac.pseudo)) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme +
    geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("GC-content of Spike In") +
    ggtitle("long Poly A pulldown \n Rep 1 Fraction Reads vs GC Content\nFrac = (rep1 + pCount)/(rep1 + rep2 + 2*pCount)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/fracError-by-gcContent.pdf"), height=12,width=5)
  
  
  # determine the relationship between fraction of reads(expected value = 1/2) and concentration
  #df.nuc.reads <- df.nuc[which(df.nuc$variable == "reads"),]
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(rep1.frac.pseudo ~ Pool14nmol.ul, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(x=Pool14nmol.ul, rep1.frac.pseudo)) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme +
    geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("Conc. of Spike In") +
    ggtitle("long Poly A pulldown \n Rep 1 Fraction Reads vs [spikeIn]\nFrac = (rep1 + pCount)/(rep1 + rep2 + 2*pCount)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/fracError-by-spikeConc.pdf"), height=12,width=5)
  
  ggplot(df.cytNuc.reads, aes(x=log10(Pool14nmol.ul), rep1.frac.pseudo, color=log10(value.rep1))) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme +
    geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("Conc. of Spike In") +
    ggtitle("long Poly A pulldown \n Rep 1 Fraction Reads vs log10[spikeIn]\nFrac = (rep1 + pCount)/(rep1 + rep2 + 2*pCount)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/fracError-by-logSpikeConc.pdf"), height=12,width=5)
  
  
  
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(log10(value.ave) ~ log10(Pool14nmol.ul), data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(x=log10(Pool14nmol.ul),color=percentGC, y=log10(value.ave))) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme + 
    geom_abline(slope=1,intercept=0,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("Spike-In nmol/ul") +
    ggtitle("long Poly A pulldown \n spike In read count vs.[spike In] ") #$#+ geom_smooth(method="lm")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/totalReads-by-spikeConc-colorByGC.pdf"), height=12,width=5)
  
  ggplot(df.cytNuc.reads, aes(x=log10(Pool14nmol.ul),color=nt, y=log10(value.ave))) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme + 
    geom_abline(slope=1,intercept=0,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("Spike-In nmol/ul") +
    ggtitle("long Poly A pulldown \n spike In read count vs.[spike In] ") #$#+ geom_smooth(method="lm")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/totalReads-by-spikeConc-colorByNT.pdf"), height=12,width=5)
  
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(value.ave ~ percentGC, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(x=percentGC, value.ave)) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme + scale_y_log10() + 
    #geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("GC Content") +
    ggtitle("long Poly A pulldown \n spike In read count vs. GC content ") 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/totalReads-by-gcContent.pdf"), height=12,width=5)
  
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- lm(value.ave ~ nt, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  ggplot(df.cytNuc.reads, aes(x=nt, value.ave)) + 
    geom_point() + #geom_smooth(method="lm")+
    facet_grid(cell~localization) + thisTheme + scale_y_log10() + 
    #geom_abline(slope=0,intercept=0.5,color="red") +
    geom_abline(data=coefs, aes(intercept=a, slope=b),color="green") + xlab("GC nucleotide len") +
    ggtitle("long Poly A pulldown \n spike In read count vs. nucleotide length ") 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/totalReads-by-nt.pdf"), height=12,width=5)
  
  # determine the relationship between fraction of reads(expected value = 1/2) and nucleotides
  df.nuc.reads <- df.nuc[which(df.nuc$variable == "reads"),]
  df.nuc.reads$localization = "nucleus"
  df.cyt.reads <- df.cyt[which(df.cyt$variable == "reads"),]
  df.cyt.reads$localization = "cytosol"
  df.reads <- rbind(df.cyt.reads,df.nuc.reads)
  df.cytNuc.reads <- df.cytNuc[which(df.cytNuc$variable == "reads"),]
  
  coefs <- ddply(df.cytNuc.reads, .(cell,localization), function(df) {
    m <- glm( rep1.frac.pseudo ~ nt, data=df)
    data.frame(a = coef(m)[1], b = coef(m)[2])
  })
  
  
}

dummy <- function(){

df.cytNuc 
df.cytNuc.rpkm <- df.cytNuc[which(df.cytNuc$variable =="transcriptTotalReads"),]
df.cytNuc.rpkm.v1_0 <- df.cytNuc.rpkm[which(df.cytNuc.rpkm$value.rep1 == 0),]
df.cytNuc.rpkm.v2_0 <- df.cytNuc.rpkm[which(df.cytNuc.rpkm$value.rep2 == 0),]
df.cytNuc.rpkm.v1_0$rep = 1
df.cytNuc.rpkm.v1_0$value = df.cytNuc.rpkm.v1_0$value.rep2
df.cytNuc.rpkm.v2_0$rep = 2
df.cytNuc.rpkm.v2_0$value = df.cytNuc.rpkm.v2_0$value.rep1
cols<- c("cell", "gene_id", "localization", "rnaExtract", "rep", "value", "variable")
df.comb <- rbind(df.cytNuc.rpkm.v1_0[cols],df.cytNuc.rpkm.v2_0[cols])

ggplot(df.comb,aes(log10(value)))+geom_bar()+
  facet_grid(cell~localization,scale="free") + thisTheme 
}

plotDifferenceBetweenRepsFlux <- function(file =getFullPath("/data/fluxCapData-lpa-proc.tab") )  {
  
  df.cytNuc1 <- read.csv(file=file, stringsAsFactors=FALSE, sep ="\t")
  
  
  df.cytNuc.rpkmSpike <- df.cytNuc1[which(df.cytNuc1$variable == "transcriptTotalRPKM_spikeIn"),]
  #df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.rpkm <- df.cytNuc1[which(df.cytNuc1$variable == "transTotalRPKM"),]
  #df.cytNuc.rpkm[which(df.cytNuc.rpkm$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.rpkm[which(df.cytNuc.rpkm$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.rpkm80 <- df.cytNuc1[which(df.cytNuc1$variable == "RPKM_80norm"),]
  #transcriptTotalConc
  df.cytNuc.conc <- df.cytNuc1[which(df.cytNuc1$variable == "transcriptTotalConc"),]
  
  
  # rep vs. rep LOG
  ggplot(df.cytNuc.conc, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nconcentration from spike In values n\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("concentration(nm/ul) rep1") + ylab("concentration(nm/ul) rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/conc-vsReps.png"), height=12,width=5)
  

  ggplot(df.cytNuc.rpkmSpike, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by Spike In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm. rep1") + ylab("RPKM spike norm. rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkmSpikeIn-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm80, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM 80 norm. rep1") + ylab("RPKM 80 norm. rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm80-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM rep1") + ylab("RPKM rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm-vsReps.png"), height=12,width=5)
  
  
  # rep vs. rep LOG
  ggplot(df.cytNuc.conc, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nconcentration from Spike In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("nmol/ul rep1") + ylab("nmol/ul rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/conc-log-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.rpkmSpike, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by Spike In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm.rep1") + ylab("RPKM spike norm.rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkmSpikeIn-log-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.rpkm80, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM 80 norm.rep1") + ylab("RPKM 80 norm rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm80-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM.rep1") + ylab("RPKM.rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm-log-vsReps.png"), height=12,width=5)
  
  
  #spike in only 
  
  ggplot(df.cytNuc.conc[which(df.cytNuc.conc$isSpikeIn == 1),], aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nconcentration from Spike In\nLongPolyA only\nSpike in only")+ 
    geom_abline(slope=1,intercept=0) + xlab("predicted nmol/ul rep1") + ylab("predicted nmol/ul rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/conc-spikeOnly-log-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$isSpikeIn == 1),], aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by Spike In\nLongPolyA only\nSpike in only")+ 
    geom_abline(slope=1,intercept=0) + xlab("RPKM spike norm.rep1") + ylab("RPKM spike norm.rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkmSpikeIn-spikeOnly-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm80[which(df.cytNuc.rpkm80$isSpikeIn == 1),], aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by sum of middle 80\nLongPolyA only\nSpike in only")+ 
    geom_abline(slope=1,intercept=0) + xlab("RPKM 80 norm.rep1") + ylab("RPKM 80 norm.rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm80-spikeOnly-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm[which(df.cytNuc.rpkm$isSpikeIn == 1),], aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM of Spike Ins\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM.rep1") + ylab("RPKM.rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm-spikeOnly-log-vsReps.png"), height=12,width=5)
  
  
  #mRNA  only 
  ggplot(df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$region == "mRNA"),], 
         aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + 
    geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by Spike In\nLongPolyA only\nmRNA only")+ 
    geom_abline(slope=1,intercept=0)+ 
    xlab("RPKM spike norm. rep1") + ylab("RPKM spike norm rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkmSpikeIn-mRNAOnly-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm[which(df.cytNuc.rpkm$region == "mRNA"),],
         aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM of Spike Ins\nmRNA only\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM.rep1") + ylab("RPKM.rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm-mRNAOnly-log-vsReps.png"), height=12,width=5)
  
  #lncRNA
  ggplot(df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$region == "lncRNA"),],
         aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM normalized by Spike In\nlncRNA only\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm. rep1") + ylab("RPKM spike norm. rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkmSpikeIn-lncOnly-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm[which(df.cytNuc.rpkm$region == "lncRNA"),], 
         aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("Flux capicitor\nRPKM of Spike Ins\nlncRNA only\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM.rep1") + ylab("RPKM.rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/rpkm-lncOnly-log-vsReps.png"), height=12,width=5)
  
  
  m.df <- as.data.frame(group_by(df.cytNuc1, cell, localization,variable,region) %.%
                          
                          summarize(sum.rep1 = sum(value.rep1),
                                    sum.rep2 = sum(value.rep2),
                                    expr.rep1 = sum(value.rep1 > 0 ),
                                    expr.rep2 = sum(value.rep2 > 0 )))
  
  colnames(m.df) <- c("cell", "localization", "measure", "region", "sum.rep1", "sum.rep2", "expr.rep1", "expr.rep2")
  melt.df <- melt(m.df, id.vars=c("cell", "localization", "measure","region"))
  m.df$frac.rep1 = with(m.df, (sum.rep1)/(sum.rep1 + sum.rep2))
  
  
  ggplot(melt.df[which(melt.df$variable %in% c("expr.rep1", "expr.rep2") & melt.df$measure == "transTotalRPKM"),], aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + 
    facet_grid(cell ~localization) + xlab("replicate") + ylab("count")+
    ggtitle("Flux capicitor\nnumber of expressed genes")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/expr-cytNuc.png"), height=12,width=5)
  
  # localization vs. cell 
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) +  
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") + ylim(0,1) +
    facet_grid(localization~cell) + 
    scale_x_discrete(limits=c("transcriptTotalReads", "transTotalRPKM","RPKM_80norm","transcriptTotalConc", "transcriptTotalRPKM_spikeIn"),
                     labels=c("reads", "RPKM","RPKM_80","[conc]", "RPKM_norm")) +
    thisTheme + ggtitle("Flux capicitor\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cytNuc-region.png"), height=5,width=12)
  
  # combined
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") + ylim(0,1) +
    
    scale_x_discrete(limits=c("transcriptTotalReads", "transTotalRPKM","RPKM_80norm","transcriptTotalConc", "transcriptTotalRPKM_spikeIn"),
                     labels=c("reads", "RPKM","RPKM_80","[conc]", "RPKM_norm")) +
    thisTheme2 + ggtitle("Flux capicitor\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cytNuc-all-combined.png"), height=5,width=10)
  
 
  ggplot(m.df, aes(x=measure,y=frac.rep1)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") + ylim(0,1) +
      scale_x_discrete(limits=c("transcriptTotalReads", "transTotalRPKM","RPKM_80norm","transcriptTotalConc", "transcriptTotalRPKM_spikeIn"),
                     labels=c("reads", "RPKM","RPKM_80","[conc]", "RPKM_norm")) +
    thisTheme2 + ggtitle("Flux capicitor\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cytNuc-all-combined-join.png"), height=5,width=10)
  
  
  #seperate by localization only
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(~localization) + ylim(0,1) +
    scale_x_discrete(limits=c("transcriptTotalReads", "transTotalRPKM","RPKM_80norm","transcriptTotalConc", "transcriptTotalRPKM_spikeIn"),
                     labels=c("reads", "RPKM","RPKM_80","[conc]" ,"RPKM_norm")) +
    thisTheme2 + ggtitle("Flux capicitor\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMnorm\nfrac.rep1=(rep1)/(rep1 + rep2)")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cytNuc-combined.png"), height=6,width=10)
  
  
  
  
}

plotBetweenReplicatesSepByLocal <- function(file=getFullPath("/data/fluxCapData-lpa-proc.tab")){  
  
  df.cytNuc <- read.csv(file=file, stringsAsFactors=FALSE, sep ="\t")
  
  m.df <- as.data.frame(group_by(df.cytNuc, cell, localization,variable,region) %.%
                          
                          summarize(sum.rep1 = sum(value.rep1),
                                    sum.rep2 = sum(value.rep2)))
  
  colnames(m.df) <- c("cell", "localization", "measure", "region", "sum.rep1", "sum.rep2")
  melt.df <- melt(m.df, id.vars=c("cell", "localization", "measure","region"))
  m.df$frac.rep1 = with(m.df, (sum.rep1)/(sum.rep1 + sum.rep2))
  
  
  #all regions
  ggplot(subset(melt.df,localization=="cytosol"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure ~ cell,scale="free_y") + 
    ggtitle("Read counts for cytosol\nAll transcripts")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cytosol.png"), height=5,width=12)
  
  ggplot(subset(melt.df,localization=="nucleus"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for nucleus\nAll transcripts")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-nucleus.png"), height=5,width=12)
  
  #mRNA
  ggplot(subset(melt.df,localization=="cytosol" & region == "mRNA"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for cytosol\nmRNA only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cyt-mRNA.png"), height=5,width=12)
  
  ggplot(subset(melt.df,localization=="nucleus" & region == "mRNA"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for nucleus\nmRNA only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-nuc-mRNA.png"), height=5,width=12)
  
  #lncRNA
  ggplot(subset(melt.df,localization=="cytosol" & region == "lncRNA"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for cytosol\nlncRNA only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cyt-lncRNA.png"), height=5,width=12)
  
  ggplot(subset(melt.df,localization=="nucleus" & region == "lncRNA"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for nucleus\nlncRNA only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-nuc-lncRNA.png"), height=5,width=12)
  
  #other
  ggplot(subset(melt.df,localization=="cytosol" & region == "other"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for cytosol\nother transcripts only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cyt-other.png"), height=5,width=12)
  
  ggplot(subset(melt.df,localization=="nucleus" & region == "other"), aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for nucleus\nother transcripts only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-nuc-other.png"), height=5,width=12)
  
  m1.df <- as.data.frame(group_by(df.cytNuc, cell, localization,variable,isSpikeIn) %.%
                           
                           summarize(sum.rep1 = sum(value.rep1),
                                     sum.rep2 = sum(value.rep2)))
  
  colnames(m1.df) <- c("cell", "localization", "measure", "region", "sum.rep1", "sum.rep2")
  melt1.df <- melt(m1.df, id.vars=c("cell", "localization", "measure","region"))
  
  ggplot(subset(melt1.df,localization=="nucleus" & region == 1), aes(x=variable,y=value)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for nucleus\nSpike-In only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-nuc-spikeIn.png"), height=5,width=12)
  
  ggplot(subset(melt1.df,localization=="cytosol" & region == 1), aes(x=variable,y=value)) + 
    geom_bar(stat="identity") + facet_grid(measure~cell,scale="free_y") + 
    ggtitle("Read counts for cytosol\nSpike-In only")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/readCount-cyt-spikeIn.png"),height=5,width=12)
}


compareRsemFlux <- function(){
  
  df.cytNuc.fpkm$prog <- "RSEM"
  df.cytNuc.rpkm$prog <- "Flux"
  d <- rbind(df.cytNuc.fpkm, df.cytNuc.rpkm)
  m <- merge(df.cytNuc.fpkm[c("gene_id", "value.ave", "cell", "localization", "rnaExtract")],
             df.cytNuc.rpkm[c("gene_id", "value.ave", "cell", "localization", "rnaExtract")],
             by =c("gene_id", "cell", "localization", "rnaExtract"),
             suffixes=c(".rsem", ".flux"))

  ggplot(m, aes(log10(value.ave.flux), log10(value.ave.rsem))) + geom_point() +
    facet_grid(cell ~ localization)
  
  
}




plotReadDistribution <- function(){
  tool <- "Flux Capacitor "
  rpkmVarName <- "transTotalRPKM"
  saveDir <- getFullPath("plots/rnaExpr/mappedReads/starSpikeIn-ERCC/cellDistro/")
  file=getFullPath("/data/fluxCapData-lpa-proc.tab")
  
  #df.cytNuc <- read.csv(file=getFullPath("/data/fluxCapData-lpa-proc.tab"),sep="\t")
  #df.cytNuc <- df.cytNuc[which(df.cytNuc$cell == "K562" & df.cytNuc$localization == "cytosol"),]
  #exportAsTable(file=getFullPath("/data/fluxCapData-K562-lpa-proc.tab"), df=df.cytNuc)
  df.cytNuc.total <- read.csv(file=file,sep="\t")
  if(!file.exists(saveDir)){dir.create(saveDir)}
  
  
  a.df <- as.data.frame(expand.grid(cell=unique(df.cytNuc.total$cell),localization=unique(df.cytNuc.total$localization)))
  
  for(cell in unique(a.df$cell)){
    for(localization in unique(a.df$localization)){
      
      df.cytNuc <- df.cytNuc.total[which(df.cytNuc.total$cell == cell & df.cytNuc.total$localization == localization),]
      df.cytNuc$exprCat <- exprOneBothNone(df=df.cytNuc, "value.rep1", "value.rep2")
      df.cytNuc.rpkm <- df.cytNuc[which(df.cytNuc$variable == rpkmVarName),]
      rep1.expressed <- dim(df.cytNuc.rpkm[which(df.cytNuc.rpkm$value.rep1 > 0),])[1]
      rep2.expressed <- dim(df.cytNuc.rpkm[which(df.cytNuc.rpkm$value.rep2 > 0),])[1]
      
      title.vec <- paste0(tool,"distribution\n",cell," ",localization, " replicate ",c(1,2),"\n",
                          "(expressed/total) = ",c(rep1.expressed,rep2.expressed)," / ",dim(df.cytNuc.rpkm)[1],"\n",
                          "\nverticle lines are percentile by 5%(5%-100%)\nblue is 20/80,green=50")
      save.file.vec <- paste0(saveDir,cell,"-rpkmDistro-Rep",
                              c(1,2),"-",localization,"-RPKM-distro")
      
      
      #ggplot(df.cyt, aes(log10(value.rep1), log10(value.rep2)))+geom_point()+facet_wrap(~variable)
      #df.cytNuc.normSum <- df.cytNuc[which(df.cytNuc$variable == "transTotalRPKM"),]
      #df.cytNuc.normSum$variable <- "RPKM_normBySum"
      #df.cytNuc.normSum$value.rep1 <- df.cytNuc.normSum$value.rep1 / sum(df.cytNuc.normSum$value.rep1)
      #df.cytNuc.normSum$value.rep2 <- df.cytNuc.normSum$value.rep2 / sum(df.cytNuc.normSum$value.rep2)
      #df.cytNuc <- rbind(df.cytNuc, df.cytNuc.normSum)
      #ggplot(df.cytNuc, aes(log10(value.rep1), log10(value.rep2)))+geom_point()+facet_wrap(~variable) + geom_abline(intercept=0,slope=1)
      
      
      ints.rep1 <- sapply(1:20/20, function(x)quantile(df.cytNuc.rpkm$value.rep1[which(df.cytNuc.rpkm$value.rep1 > 0)], x))
      rep1.min <- min(df.cytNuc.rpkm$value.rep1[which(df.cytNuc.rpkm$value.rep1>0)])
      ggplot(df.cytNuc.rpkm, aes(log10(value.rep1),fill=exprCat)) + geom_bar() +
        geom_vline(xintercept =log10(c(rep1.min,ints.rep1)))+
        geom_vline(xintercept =log10(ints.rep1[10]),color="green")+
        geom_vline(xintercept =log10(ints.rep1[c(4,16)]),color="blue") +xlab("log10 RPKM")+
        ggtitle(title.vec[1]) + theme_bw()
      ggsave(paste0(save.file.vec[1],".png"), height=7,width=7)
      
      df.cytNuc.rpkm$value.rep1 <- ifelse( df.cytNuc.rpkm$value.rep1 == 0, rep1.min,df.cytNuc.rpkm$value.rep1 )
      ggplot(df.cytNuc.rpkm, aes(log10(value.rep1),fill=exprCat)) + geom_bar() +
        geom_vline(xintercept =log10(c(rep1.min,ints.rep1)))+
        geom_vline(xintercept =log10(ints.rep1[10]),color="green")+
        geom_vline(xintercept =log10(ints.rep1[c(4,16)]),color="blue") +xlab("log10 RPKM")+
        ggtitle(paste(title.vec[1],"\nRPKM=0 displayed @ minimum")) + theme_bw()
      ggsave(paste0(save.file.vec[1],"minDisplayed.png"), height=7,width=7)
      
      
      
      ints.rep2 <- sapply(1:20/20, function(x)quantile(df.cytNuc.rpkm$value.rep2[which(df.cytNuc.rpkm$value.rep2 > 0)], x))  
      rep2.min <- min(df.cytNuc.rpkm$value.rep2[which(df.cytNuc.rpkm$value.rep2>0)])
      
      ggplot(df.cytNuc.rpkm, aes(log10(value.rep2),fill=exprCat)) + geom_bar() +
        geom_vline(xintercept =log10(c(rep2.min,ints.rep2)))+
        geom_vline(xintercept =log10(ints.rep2[10]),color="green")+
        geom_vline(xintercept =log10(ints.rep2[c(4,16)]),color="blue") +xlab("log10 RPKM")+
        ggtitle(title.vec[2])+ theme_bw()
      ggsave(paste0(save.file.vec[2],".png"), height=7,width=7)
      
      
      df.cytNuc.rpkm$value.rep2 <- ifelse( df.cytNuc.rpkm$value.rep2 == 0, rep2.min,df.cytNuc.rpkm$value.rep2 )
      ggplot(df.cytNuc.rpkm, aes(log10(value.rep2),fill=exprCat)) + geom_bar() +
        geom_vline(xintercept =log10(c(rep2.min,ints.rep2)))+
        geom_vline(xintercept =log10(ints.rep2[10]),color="green")+
        geom_vline(xintercept =log10(ints.rep2[c(4,16)]),color="blue") +xlab("log10 RPKM")+
        ggtitle(paste(title.vec[2],"\nRPKM=0 displayed @ minimum")) + theme_bw()
      ggsave(paste0(save.file.vec[2],"minDisplayed.png"), height=7,width=7)
      
    }
  }
  
}






