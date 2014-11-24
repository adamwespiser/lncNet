


doitVersion2 <- function(){
  v2.infile = getFullPath("/data/rsemCapData-v2-lpa-proc.tab")
  v2.outdir = getFullPath("plots/rnaExpr/mappedReads/RSEM-version2/")
  plotDifferenceBetweenRepsRSEM(infile=v2.infile,outdir=v2.outdir)
  plotRSEMcytFrac(infile=v2.infile,outdir=v2.outdir)
  
  plotRSEMreadDistro(inFile=getFullPath("/data/rsemCapData-v2-readDistro.tab"),
                    outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM-version2/readDistro/"),
                    summaryPlotTitle=" bowtie -> RSEM\nmapped reads")
  
}



doitRSEM_RSEM <- function(){
  v2.infile =getFullPath("/data/rsemRsemmapped-v1-lpa-proc.tab")
  v2.outdir = getFullPath("plots/rnaExpr/mappedReads/RSEM-rsemMapped/")
 
  plotDifferenceBetweenRepsRSEM(infile=v2.infile,outdir=v2.outdir,plotMsg="RSEM -> ")
  
  plotRSEMcytFrac(infile=v2.infile,outdir=v2.outdir)
  
  plotRSEMreadDistro(inFile= getFullPath("/data/rsemRsemmapped-v1-readDistro.tab"),
                     outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM-rsemMapped/readDistro/"),
                     summaryPlotTitle=" RSEM(bowtie) -> RSEM\nmapped reads")
  
}


doitSTAR_RSEM <- function(){
  v2.infile =getFullPath("/data/rsemSTARmapped-v1-lpa-proc.tab")
  v2.outdir = getFullPath("plots/rnaExpr/mappedReads/RSEM-STAR-mapped/")
  plotDifferenceBetweenRepsRSEM(infile=v2.infile,outdir=v2.outdir,plotMsg="STAR -> ")
  plotRSEMcytFrac(infile=v2.infile,outdir=v2.outdir)
  
  plotRSEMreadDistro(inFile=getFullPath("/data/rsemSTARmapped-v1-readDistro.tab"),
                     outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM-STAR-mapped/readDistro/"),
                     summaryPlotTitle=" STAR -> RSEM\nmapped reads")
  
}

doitBowtie_RSEM <- function(){
  v2.infile =getFullPath("/data/rsemBowtiemapped-v1-lpa-proc.tab")
  v2.outdir = getFullPath("plots/rnaExpr/mappedReads/RSEM-bowtie-mapped/")
  plotDifferenceBetweenRepsRSEM(infile=v2.infile,outdir=v2.outdir,plotMsg="bowtie -> ")
  plotRSEMcytFrac(infile=v2.infile,outdir=v2.outdir)
  
  plotRSEMreadDistro(inFile=getFullPath("/data/rsemBowtiemapped-v1-readDistro.tab"),
                     outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM-bowtie-mapped/readDistro/"),
                     summaryPlotTitle=" STAR -> RSEM\nmapped reads")
  
}



plotDifferenceBetweenRepsRSEM <- function(infile=getFullPath("/data/rsemCapData-lpa-proc.tab"),
                                          outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM/"),
                                          plotMsg="")  {
  
  stopifnot(file.exists(infile))
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  
  df.cytNuc <- read.csv(file=infile, stringsAsFactors=FALSE, sep ="\t")
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
    ggtitle(paste0(plotMsg,"RSEM\nTPM \nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("TPM. rep1") + ylab("TPM rep2")
  ggsave(paste0(outdir,"tpm-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.fpkm80, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle(paste0(plotMsg,"RSEM\nFPKM normalized by sum of middle 80 In\nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM 80 norm. rep1") + ylab("FPKM 80 norm. rep2")
  ggsave(paste0(outdir,"fpkm80-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(x=value.rep1,y=value.rep2,color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle(paste0(plotMsg,"RSEM\nFPKM\nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("FPKM rep1") + ylab("FPKM rep2")
  ggsave(paste0(outdir,"fpkm-vsReps.png"), height=12,width=5)
  
  
  # rep vs. rep LOG
  ggplot(df.cytNuc.tpm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle(paste0(plotMsg,"RSEM\nTPM n\nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("log10 of TPM rep1") + ylab("log10 of TPM rep2")
  ggsave(paste0(outdir,"tpm-log-vsReps.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.fpkm80, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle(paste0(plotMsg,"RSEM\nFPKM normalized by sum of middle 80 In\nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("log10 of  FPKM 80 norm. rep1") + ylab("log10 of FPKM 80 norm. rep2")
  ggsave(paste0(outdir,"fpkm80-log-vsReps.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
    facet_grid(cell~localization,scale="free") + thisTheme +
    ggtitle(paste0(plotMsg,"RSEM\nFPKM \nLongPolyA only"))+ 
    geom_abline(slope=1,intercept=0)+ xlab("log10 of FPKM rep1") + ylab("log10 of FPKM rep2")
  ggsave(paste0(outdir,"fpkm-log-vsReps.png"), height=12,width=5)
  
  
  # Plot cutFraction
  
  
  
  
  
  #mRNA  only 
  #   ggplot(df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$region == "mRNA"),], 
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + 
  #     geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nFPKM normalized by Spike In\nLongPolyA only\nmRNA only")+ 
  #     geom_abline(slope=1,intercept=0)+ 
  #     xlab("RPKM spike norm.") + ylab("RPKM spike norm")
  #   ggsave(paste0(outdir,"rpkmSpikeIn-mRNAOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   ggplot(df.cytNuc.fpkm[which(df.cytNuc.fpkm$region == "mRNA"),],
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM of Spike Ins\nmRNA only\nLongPolyA only")+ 
  #     geom_abline(slope=1,intercept=0)+ xlab("RPKM") + ylab("RPKM ")
  #   ggsave(paste0(outdir,"rpkm-mRNAOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   #lncRNA
  #   ggplot(df.cytNuc.fpkmSpike[which(df.cytNuc.fpkmSpike$region == "lncRNA"),],
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM normalized by Spike In\nlncRNA only\nLongPolyA only")+ 
  #     geom_abline(slope=1,intercept=0)+ xlab("RPKM spike norm.") + ylab("RPKM spike norm")
  #   ggsave(paste0(outdir,"rpkmSpikeIn-lncOnly-log-vsReps.png"), height=12,width=5)
  #   
  #   ggplot(df.cytNuc.fpkm[which(df.cytNuc.fpkm$region == "lncRNA"),], 
  #          aes(x=log10(value.rep1),y=log10(value.rep2),color=region)) + geom_point() + 
  #     facet_grid(cell~localization,scale="free") + thisTheme +
  #     ggtitle("RSEM\nRPKM of Spike Ins\nlncRNA only\nLongPolyA only")+ 
  #     geom_abline(slope=1,intercept=0)+ xlab("RPKM") + ylab("RPKM ")
  #   ggsave(paste0(outdir,"rpkm-lncOnly-log-vsReps.png"), height=12,width=5)
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
    geom_bar(stat="identity") +  xlab("replicates") + ylab("count") +
    facet_grid(cell ~localization) +
  ggtitle(paste0(plotMsg,"RSEM\nnumber of genes with mapped reads"))
  ggsave(paste0(outdir,"expr-cytNuc.png"), height=12,width=5)
  # localization vs. cell 
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) +  
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(localization~cell) + 
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme + ggtitle(paste0(plotMsg,"RSEM\nfraction of cytosol & nucleus \nFPKM/TPM/FPKM_80\nfrac.rep1=(rep1)/(rep1 + rep2)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"readCount-cytNuc-region.png"), height=5,width=12)
  
  # combined
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle(paste0(plotMsg,"RSEM\nfraction of cytosol & nucleus \nFPKM/TPM/FPKM_80\n frac.rep1=(rep1)/(rep1 + rep2)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"readCount-cytNuc-all-combined.png"), height=5,width=10)
  
  ggplot(m.df, aes(x=measure,y=frac.rep1)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle(paste0(plotMsg,"RSEM\nfraction of cytosol & nucleus \nFPKM/TPM/FPKM_80\nfrac.rep1=(rep1)/(rep1 + rep2)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"readCount-cytNuc-all-combined-join.png"), height=5,width=10)
  
  
  #seperate by localization only
  ggplot(m.df, aes(x=measure,y=frac.rep1,color=region)) + 
    geom_boxplot() + geom_abline(slope=0,intercept=1/2,color="red") +
    facet_grid(~localization) + 
    scale_x_discrete(limits=c("FPKM","TPM","FPKM_80norm"),
                     labels=c("FPKM", "TPM","FPKM_80")) + ylim(0,1) +
    thisTheme2 + ggtitle(paste0(plotMsg,"RSEM\nfraction of cytosol & nucleus \nFPKM/TPM/FPKM_80\n frac.rep1=(rep1)/(rep1 + rep2)"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(paste0(outdir,"readCount-cytNuc-combined.png"), height=6,width=10)
  
}

plotReadDistributionRSEM <- function(outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM/"),
                                     infile=getFullPath("/data/rsemCapData-lpa-proc.tab"),
                                     tool="RSEM"){
  if(!file.exists(outdir)){dir.create(path=outdir,recursive=TRUE)}
  stopifnot(file.exists(infile))
  
  
  
  rpkmVarName <- "FPKM"
  saveDir <- paste0(outdir,"cellDistro/")
  
  file = infile
  #df.cytNuc <- read.csv(file=getFullPath("/data/fluxCapData-lpa-proc.tab"),sep="\t")
  #df.cytNuc <- df.cytNuc[which(df.cytNuc$cell == "K562" & df.cytNuc$localization == "cytosol"),]
  #exportAsTable(file=getFullPath("/data/fluxCapData-K562-lpa-proc.tab"), df=df.cytNuc)
  df.cytNuc.total <- read.csv(file=file,sep="\t")
  #if(!file.exists(saveDir)){dir.create(saveDir)}
  
  
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

plotRSEMcytFrac <- function(infile=getFullPath("/data/rsemCapData-lpa-proc.tab"),
                            outdir=getFullPath("/plots/rnaExpr/mappedReads/RSEM/")){
  
  if(!file.exists(outdir)){dir.create(outdir,recursive=TRUE)}
  stopifnot(file.exists(infile))
  
  df.cytNuc <- read.csv(file=infile, stringsAsFactors=FALSE, sep ="\t")
  
  
  df.cyt <- df.cytNuc[which(df.cytNuc$localization == "cytosol"),]
  df.nuc <- df.cytNuc[which(df.cytNuc$localization == "nucleus"),]
  df.cytNuc1 <- merge(df.cyt,df.nuc,by=c("gene_id","cell","variable"),suffixes=c(".cyt",".nuc"))
  df.cytNuc1$cytFracPseudo <- with(df.cytNuc1, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  df.cytNuc1$cytFrac <- with(df.cytNuc1, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  
  
  df.cytNuc.fpkm <- df.cytNuc1[which(df.cytNuc1$variable =="FPKM"),]
  #df.cytNuc.fpkm[which(df.cytNuc.fpkm$gene_id %in% pc),"region"] <- "mRNA"
  #df.cytNuc.fpkm[which(df.cytNuc.fpkm$gene_id %in% lnc),"region"] <- "lncRNA"
  df.cytNuc.fpkm80<- df.cytNuc1[which(df.cytNuc1$variable =="FPKM_80norm"),]
  #transcriptTotalConc
  df.cytNuc.tpm<- df.cytNuc1[which(df.cytNuc1$variable =="TPM"),]
  
  
  ggplot(df.cytNuc.tpm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nTPM: cyt/(nuc + cyt)")
  ggsave(paste0(outdir,"cytFrac-tpm.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nFPKM: cyt/(nuc + cyt)")
  ggsave(paste0(outdir,"cytFrac-fpkm.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.fpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFrac,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nFPKM80): cyt/(nuc + cyt)")
  ggsave(paste0(outdir,"cytFrac-fpkm80.png"), height=12,width=5)
  
  
  #PLOT pseudo cytFrac
  ggplot(df.cytNuc.tpm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nTPM: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(paste0(outdir,"cytFracPseudo-tpm.png"), height=12,width=5)
  
  ggplot(df.cytNuc.fpkm, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nFPKM: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(paste0(outdir,"cytFracPseudo-fpkm.png"), height=12,width=5)
  
  
  ggplot(df.cytNuc.fpkm80, aes(y=log10(value.ave.cyt*2 + value.ave.nuc*2),x=cytFracPseudo,color=factor(region.cyt)))+
    geom_density2d() + theme_bw() + thisTheme + 
    facet_grid(cell~region.cyt)+
    ggtitle("RSEM:  \nFraction of Cytosolic RNA-seq expr\nFPKM80: cytPseudo/(nucPseudo + cytPseudo)")
  ggsave(paste0(outdir,"cytFracPseudo-fpkm80.png"), height=12,width=5)
  
}

plotRSEMreadDistro <- function(inFile=getFullPath("/data/rsemCapData-v2-readDistro.tab"),
                               outdir=getFullPath("plots/rnaExpr/mappedReads/RSEM-version2/readDistro/"),
                               summaryPlotTitle=" RSEM: mapped reads ",
                               removeCells = c("GM12878", "HeLa-S3", "HUVEC", "MCF-7") ){
  
  if(!file.exists(outdir)){dir.create(outdir,recursive=TRUE)}
  rdCnt <- read.csv(file=inFile, stringsAsFactors=FALSE, sep ="\t")
  rdCnt$rep <- ifelse(rdCnt$replicate > 2, rdCnt$replicate - 2, rdCnt$replicate)
  rdCnt$readCount <- rdCnt$count
  
  if(length(removeCells) > 0){
    rdCnt <- rdCnt[-which(rdCnt$cell %in% removeCells),]
  }
  
  if(FALSE){
  for(rna in c("longPolyA")){
    for(loc in c("cytosol","nucleus")){
      base = paste0(outdir,rna,"-",loc,"-")
      localTitle = paste(summaryPlotTitle,rna,loc,"\nFacet by Replicate,Celltype\nmultiplicity=genome aligns per read")
     
      ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
             aes(x=multiplicity,y=count)) + geom_line() + geom_point()+ 
        facet_grid(cell~rep) + theme_bw()+
        ggtitle(localTitle)+
      ggsave(paste0(base,"readDistro.pdf"), height=12,width=5)
      
      ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
             aes(x=multiplicity,y=count)) + geom_line() + geom_point()+ 
        facet_grid(cell~rep) + theme_bw()+
      ggtitle(localTitle) + xlim(0,20)
      ggsave(paste0(base,"readDistro-zoom.pdf"), height=12,width=5)
      
      
      ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
             aes(x=log10(multiplicity),y=count)) + geom_line() + geom_point()+ 
        facet_grid(cell~rep) + theme_bw() +
      ggtitle(localTitle) 
      ggsave(paste0(base,"readDistro-logMult.pdf"), height=12,width=5)
      
      ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
             aes(x=log10(multiplicity),y=log10(count))) + geom_line() + geom_point()+ 
        facet_grid(cell~rep) + theme_bw() +
      ggtitle(localTitle) 
      ggsave(paste0(base,"readDistro-logMultLogCount.pdf"), height=12,width=5)
      
      ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
             aes(x=multiplicity,y=log10(count))) + geom_line() + geom_point()+ 
        facet_grid(cell~rep) + theme_bw() +
        ggtitle(localTitle) 
      ggsave(paste0(base,"readDistro-LogCount.pdf"), height=12,width=5)
    }  
    
  }
  }
  byExp <- as.data.frame(group_by(rdCnt, rep,Exp,replicate,rnaExtract,localization,cell) %>% summarise(allReads = sum(count)))
  byExp.cyt <- subset(byExp, localization=="cytosol")
  byExp.nuc <- subset(byExp, localization=="nucleus")
 
  ggplot(byExp, aes(x=cell,y=allReads,fill=as.factor(rep))) + 
    geom_bar(stat="identity",position="dodge") +
    xlab("RNA sequencing Expr.") + ylab("Reads") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1) )+
    facet_grid(localization ~ . ) +
    ggtitle(paste0(summaryPlotTitle,"\nReads used by RSEM\nmapped by star"))
  ggsave(paste0(outdir,"reads--Count-Used-By-RSEM.pdf"), height=5,width=8)
  
  logymax = log10(max(byExp$allReads))
  
  ggplot(byExp, aes(x=cell,y=allReads,fill=as.factor(rep))) + 
    geom_bar(stat="identity",position="dodge") +
    xlab("RNA sequencing Expr.") + ylab("log10(Reads)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1) )+
    facet_grid(localization ~ . ) +
    scale_y_log10(breaks=10^seq(from=0,to=ceiling(log10(max(byExp$allReads)))))+
    ggtitle(paste0(summaryPlotTitle,"\nReads used by RSEM\nmapped by star"))
  ggsave(paste0(outdir,"reads--log10ofCount-Used-By-RSEM.pdf"), height=5,width=8)
  
  
  
  
  
  summ <- ddply(rdCnt, .(Exp),summarise,notMapped=mean(notMapped),mapped=mean(mapped),rnaExtract=rnaExtract[1])
  ggplot(melt(summ,id.var="Exp"), aes(x=Exp,y=value,fill=variable)) + 
    geom_bar(stat="identity")+
    xlab("RNA sequencing Expr.") + ylab("Reads Mapped")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(summaryPlotTitle)
  ggsave(paste0(outdir,"readMap-Count.pdf"), height=5,width=12)
  
  ggplot(summ, aes(x=Exp,y=mapped)) + 
    geom_bar(stat="identity")+
    xlab("RNA sequencing Expr.") + ylab("Reads Mapped")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(summaryPlotTitle)+
    geom_abline(slope=0,intercept=20*(10^6))
  ggsave(paste0(outdir,"readMap-mappedOnly.pdf"), height=5,width=12)
  
  ggplot(summ, aes(x=Exp,y=mapped)) + 
    geom_bar(stat="identity")+
    xlab("RNA sequencing Expr.") + ylab("Reads Mapped")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(summaryPlotTitle)+
    geom_abline(slope=0,intercept=20*(10^6))+
    facet_grid(rnaExtract~.)
  ggsave(paste0(outdir,"readMap-mappedOnly-splitByPolyA.pdf"), height=6,width=12)
  
  ggplot(summ, aes(x=Exp,y=100*mapped/(mapped +notMapped))) + 
    geom_bar(stat="identity")+
    ggtitle(summaryPlotTitle)+
    xlab("RNA sequencing Expr.") + ylab("% Reads mapped")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylim(0,100) +
    facet_grid(rnaExtract~.)
  ggsave(paste0(outdir,"readMap-Percent-splitByPolyA.pdf"), height=6,width=12)
  
  
  
  
  ggplot(rdCnt[which(rdCnt$rnaExtract == rna & rdCnt$localization == loc ),], 
         aes(x=multiplicity,y=count)) + geom_line() + geom_point()+ 
    facet_grid(cell~rep) + theme_bw()
  
}


