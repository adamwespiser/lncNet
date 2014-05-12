

plotDifferenceBetweenRepsRPKMfromBam <- function(file=getFullPath("/data/rpkmFromBamCapData-lpa-proc.tab") )  {
  
  df.cytNuc <- read.csv(file=file, stringsAsFactors=FALSE, sep ="\t")
  #df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.rpkm <- df.cytNuc[which(df.cytNuc$variable =="RPKM"),]
  #df.cytNuc.rpkm[which(df.cytNuc.rpkm$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.rpkm[which(df.cytNuc.rpkm$gene_id %in% lnc),"region"] <- "lncRNA"
  
  df.cytNuc.rpkm80<- df.cytNuc[which(df.cytNuc$variable =="RPKM_80norm"),]
  df.cytNuc.concBySpikeIn<- df.cytNuc[which(df.cytNuc$variable =="concBySpikeIn"),]
  df.cytNuc.spikeIn_norm<- df.cytNuc[which(df.cytNuc$variable =="spikeIn_norm"),]
  
  #spikeIn_norm
  
  #transcriptTotalConc
  #df.cytNuc.tpm<- df.cytNuc[which(df.cytNuc$variable =="TPM"),]
  
  
  
  
  # rep vs. rep LOG
#   ggplot(df.cytNuc.tpm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RPKMfromBam\nTPM \nLongPolyA only")+ 
#     geom_abline(slope=1,intercept=0)+ xlab("TPM") + ylab("TPM")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/tpm-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.rpkm80, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM 80 norm. rep1") + ylab("RPKM 80 norm. rep1")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkm80-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM rep1") + ylab("RPKM rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkm-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.concBySpikeIn, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nconcentration from spikeIn glm\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("rep1:[nmol/ul]") + ylab("rep2:[nmol/ul]")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/conc-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.spikeIn_norm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM/(sum spikeIn RPKM)\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM/(sum spikeIn RPKM ) rep1") + ylab("RPKM/(sum spike In RPKM) rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkmSpikeIn-vsReps.png"), height=12,width=5)
  
  
  
  # rep vs. rep LOG
#   ggplot(df.cytNuc.tpm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
#     facet_grid(cell~localization,scale="free") + thisTheme +
#     ggtitle("RPKMfromBam\nTPM n\nLongPolyA only")+ 
#     geom_abline(slope=1,intercept=0)+ xlab("TPM") + ylab("TPM")
#   ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/tpm-log-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.rpkm80, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM normalized by sum of middle 80 In\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM 80 norm. rep1") + ylab("RPKM 80 norm. rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkm80-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.rpkm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM \nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM rep1") + ylab("RPKM rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkm-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.concBySpikeIn, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nconcentration from spike in glm\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("rep1:[nmol/ul]") + ylab("rep2:[nmol/ul]")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/conc-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.spikeIn_norm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle("RPKMfromBam\nRPKM/sum spikeIn RPKM\nLongPolyA only")+ 
    geom_abline(slope=1,intercept=0)+ xlab("RPKM/(sum spikeIn RPKM ) rep1") + ylab("RPKM/(sum spike In RPKM) rep2")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/rpkmSpikeIn-log-vsReps.png"), height=12,width=5)
  
  
  #mRNA  only 
  #   ggplot(df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$region == "mRNA"),], 
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + 
  #     geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM normalized by Spike In\nLongPolyA only\nmRNA only")+ 
  #     geom_abline(slope=1,intercept=0)+ 
  #     xlab("RPKM spike norm.") + ylab("RPKM spike norm")
  #   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkmSpikeIn-mRNAOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   ggplot(df.cytNuc.rpkm[which(df.cytNuc.rpkm$region == "mRNA"),],
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM of Spike Ins\nmRNA only\nLongPolyA only")+ 
  #     geom_abline(slope=1,intercept=0)+ xlab("RPKM") + ylab("RPKM ")
  #   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkm-mRNAOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   #lncRNA
  #   ggplot(df.cytNuc.rpkmSpike[which(df.cytNuc.rpkmSpike$region == "lncRNA"),],
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM normalized by Spike In\nlncRNA only\nLongPolyA only")+ 
  #     geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm.") + ylab("RPKM spike norm")
  #   ggsave(getFullPath("plots/rnaExpr/mappedReads/RSEM/rpkmSpikeIn-lncOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   ggplot(df.cytNuc.rpkm[which(df.cytNuc.rpkm$region == "lncRNA"),], 
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
  
  
  ggplot(melt.df[which(melt.df$variable %in% c("expr.rep1", "expr.rep2") & melt.df$measure == "RPKM"),], aes(x=variable,y=value,fill=region)) + 
    geom_bar(stat="identity") + 
    facet_grid(cell ~localization) + xlab("replicates") + ylab("count") +
    ggtitle("RPKMfromBam\nnumber of expressed genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/expr-cytNuc.png"), height=12,width=5)
  # localization vs. cell 
  # "RPKM_80norm","RPKM","reads","concBySpikeIn","spikeIn_norm
  
  
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) +  
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(localization~cell) + 
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","spikeIn_norm","concBySpikeIn"),
                     labels=c("reads","RPKM"  ,"RPKM_80","RPKM/spikeIn","Conc.")) + ylim(0,1) + 
    thisTheme + 
    ggtitle("RPKMfromBam\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMspikeIn/Conc\nfrac.rep1=(rep1)/(rep1 + rep2)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/readCount-cytNuc-region.png"), height=5,width=12)
  
  # combined
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","spikeIn_norm","concBySpikeIn"),
                     labels=c("reads","RPKM"  ,"RPKM_80","RPKM/spikeIn","Conc."))+ ylim(0,1) +
    thisTheme2 + 
    ggtitle("RPKMfromBam\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80/RPKMspikeIn/Conc\n frac.rep1=(rep1)/(rep1 + rep2)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/readCount-cytNuc-all-combined.png"), height=5,width=10)
  
  ggplot(m.df, aes(x=measure,y=frac.rep1)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","spikeIn_norm","concBySpikeIn"),
                     labels=c("reads","RPKM"  ,"RPKM_80","RPKM/spikeIn","Conc."))+ ylim(0,1) +
    thisTheme2 + 
    ggtitle("RPKMfromBam\nfraction of cytosol & nucleus\nreads/RPKM/RPKM80/RPKMspikeIn/Conc\nfrac.rep1=(rep1)/(rep1 + rep2)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/readCount-cytNuc-all-combined-join.png"), height=5,width=10)
  
  
  #seperate by localization only
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(~localization) +
    scale_x_discrete(limits=c("reads","RPKM", "RPKM_80norm","spikeIn_norm","concBySpikeIn"),
                     labels=c("reads","RPKM"  ,"RPKM_80","RPKM/spikeIn","Conc.")) + ylim(0,1) +
    thisTheme2 + ggtitle("RPKMfromBam\nfraction of cytosol & nucleus \nreads/RPKM/RPKM80\nfrac.rep1=(rep1)/(rep1 + rep2)")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/readCount-cytNuc-combined.png"), height=6,width=10)
}


plotReadDistributionRPKMfromBAM <- function(){
  tool <- "RPKMfromBam "
  rpkmVarName <- "RPKM"
  saveDir <- getFullPath("plots/rnaExpr/mappedReads/RPKMfromBam/cellDistro/")
  file=getFullPath("/data/rpkmFromBamCapData-lpa-proc.tab")  
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






