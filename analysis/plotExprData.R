homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }

source(getFullPath("analysis/gencodeInput.R"))


local.datadir <<- "/home/wespisea/data/"
local.plotdir.encodeExpr <<- getFullPath("plots/rnaExpr/ENCODE-cytVnuc/")
getEPlotPath <- function(subpath){ file.path(local.plotdir.encodeExpr, subpath) }


plotENCODErnaseq_annot <- function(){
  annot <- getAnnotFile()
  ggplot(annot, aes(cell, localization)) + geom_tile() +
    ggtitle("RNA-seq V10 mapped") + coord_flip() + theme_bw() + 
    facet_grid(. ~pulldown)
  ggsave(getEPlotPath("RNAseqDataSource.pdf"),height=7,width=7)
}


plotENCODErnaseq <- function(){
  if(!file.exists(local.plotdir.encodeExpr)){
    dir.create(local.plotdir.encodeExpr,recursive=T)
  }
  df <- cytNucRatios()
  
  ### Get number of found genes/celltype
  df$Count <- 1
  dgc.table <- tbl_df(df)
  dgc.table.cell <- group_by(dgc.table, cell,gene_type)
  df.cell.counts <- as.data.frame(summarise(dgc.table.cell,
                              num.types = n()))
  

  
  ggplot(df.cell.counts, aes(x=cell, y = num.types, fill=gene_type))+ geom_bar(stat="identity") +
    coord_flip() + theme_bw() + 
    ggtitle("genes w/ expr > 0 in (cyt,nuc)")
  ggsave(getEPlotPath("cellTypeCountByGeneType.pdf"),height=7,width=7) #1
  
  
  lncwFnc <- filter(df,biotype == "lnc") %.%
    group_by(gene_type,cell,funcLnc) %.%
    summarise(num.types = n()) 
  
  ggplot(as.data.frame(lncwFnc), aes(x=cell, y = num.types, fill=gene_type))+ geom_bar(stat="identity") +
    theme_bw() + facet_grid(funcLnc ~ . ,scale="free")+
    ggtitle("lncRNA w/ expr > 0 in (cyt,nuc)\n{0,1} == {func,unlabelled}") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(getEPlotPath("cellTypeCountByLncFunc.pdf"),height=7,width=7) #1
  
  
  
  
  funcLnc.biotype <- unique(filter(df, biotype == "lnc", funcLnc == 1) %.%
                              select(gene_id,gene_type)) %.%
    group_by(gene_type) %.%
    summarise(num.types = n())
  
  ggplot(funcLnc.biotype, aes(x=gene_type, y = num.types)) + geom_bar(stat="identity") +
    coord_flip() + theme_bw() +
    ggtitle("functional Lnc found by biotype")
  ggsave(getEPlotPath("funcLncCountByGeneType.pdf"),height=7,width=7) # 2
  
  ggplot(select(df,funcLnc==1), aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("functional Lnc found by biotype") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")
  ggsave(getEPlotPath("cytVsNuc-lncFunc.pdf"),height=7,width=7) # 3
  
  ggplot(subset(df,biotype == "lnc"), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc)))+
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("All Lnc: cyt vs. nuc\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")+
    scale_size_manual(values = c(3,4))
  ggsave(getEPlotPath("cytVsNuc-All-lnc.pdf"),height=7,width=7) # 4
  
  ggplot(subset(df,biotype == "lnc"& (RPKMsum.cyt < 500) & (RPKMsum.nuc < 500)), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc)))+
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("All Lnc: cyt vs. nuc\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")+
    scale_size_manual(values = c(3,4))
  ggsave(getEPlotPath("cytVsNuc-All-lnc-lt500.pdf"),height=7,width=7) 
  
  
  
  ggplot(subset(df,biotype=="lnc"), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~gene_type, nrow=1,scale="free") +
    ggtitle("biotype=Lnc by gene_type\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,4))
  ggsave(getEPlotPath("cytVsNuc-lncBy.pdf"),height=5,width=10) # 5
  
  ggplot(arrange(subset(df,biotype=="lnc" & RPKMsum.nuc < 500 & RPKMsum.cyt < 500),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~gene_type, nrow=1,scale="free") +
    ggtitle("biotype=Lnc by gene_type\ny=x abline,RPKMsum... < 500") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3.5))
  ggsave(getEPlotPath("cytVsNuc-lncBy-lt500.pdf"),height=5,width=10) 
  
  ggplot(df,aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~biotype, nrow=1,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-biotypeCompare.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.nuc < 500 & RPKMsum.cyt < 500),aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point(alpha=I(0.5)) +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~biotype, nrow=1,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-biotypeCompare-lt500.pdf"),height=5,width=10)
  
  
  ggplot(df,aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~biotype,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare.pdf"),height=24,width=7)
  
  ## cyt/nuc sum comparison
  ggplot(subset(df, biotype=="lnc"), aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
  
  ggplot(arrange(subset(df, biotype=="lnc" & (RPKMsum.cyt < 500) & (RPKMsum.nuc < 500)),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3))
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA-lt500.pdf"),height=24,width=7)
  
  ggplot(arrange(subset(df, biotype=="lnc" & (RPKMsum.cyt < 100) & (RPKMsum.nuc < 100)),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3)) +xlim(0,100) + ylim(0,100)
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA-lt100.pdf"),height=24,width=7)
  
  
  ## Correlation between replicates
  cyt.rsq <- as.data.frame(group_by(df,cell) %.% summarise(count = n()))
  cyt.rsq$count <- NULL
  cyt.rsq$gene_type <- "total"
  cyt.rsq$rsquaredCyt <- as.numeric(group_by(df,cell) %.% do(function(x)summary(lm(RPKM1.cyt ~ RPKM2.cyt, data=x))$r.squared))
  cyt.rsq$rsquaredNuc <- as.numeric(group_by(df,cell) %.% do(function(x)summary(lm(RPKM1.nuc ~ RPKM2.nuc, data=x))$r.squared))
 cyt.rsq.lnc <- as.data.frame(group_by(df,cell,gene_type) %.% summarise(count = n()))
  cyt.rsq.lnc$count <- NULL
  cyt.rsq.lnc$rsquaredCyt <- as.numeric(group_by(df,cell,gene_type) %.% do(function(x)summary(lm(RPKM1.cyt ~ RPKM2.cyt, data=x))$r.squared))
  cyt.rsq.lnc$rsquaredNuc <- as.numeric(group_by(df,cell,gene_type) %.% do(function(x)summary(lm(RPKM1.nuc ~ RPKM2.nuc, data=x))$r.squared))
  cytNuc.comb.rsq <- rbind(cyt.rsq.lnc,cyt.rsq)
  
  ggplot(melt(cytNuc.comb.rsq, id.var=c("cell","gene_type")), aes(x=cell,y=value,fill=value)) + geom_bar(stat="identity") +
    facet_grid(variable~gene_type,scale="free") + theme_bw() +
    xlab("R-squared for subset/total in cell type") + ylab("cell type") + 
    ggtitle("lm on RPKM1 & RPKM2 of cytosol\nlncRNA")  +coord_flip() +
    scale_fill_gradient2(low="red",mid="grey",high="green",midpoint=0.75)
  ggsave(getEPlotPath("cyt-Stats-allgenes.pdf"),height=7,width=7)

  
  
  ## cyt sum comparison
  

  ggplot(subset(df, biotype=="lnc"), 
         aes(x = RPKM1.cyt,y=RPKM2.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cytosols replicates \n lncRNA facet by genetype\ny=x abline") +
    xlab("RPKM 1 cytosol") + ylab("RPKM 2 cytosol") 
  ggsave(getEPlotPath("cyt-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
 
  ggplot(subset(df, biotype=="lnc" & (RPKM1.cyt < 500) & (RPKM2.cyt < 500)), 
         aes(x = RPKM1.cyt,y=RPKM2.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cytosols replicates \n lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("RPKM 1 cytosol") + ylab("RPKM 2 cytosol") 
  ggsave(getEPlotPath("cyt-cellVbiotypeCompare-lncRNA-lncRNA-lt500.pdf.pdf"),height=24,width=7)
  

  ## nuc sum comparison
  ggplot(subset(df, biotype=="lnc"), 
         aes(x = RPKM1.nuc,y=RPKM2.nuc)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("nucleus replicates \n lncRNA facet by genetype\ny=x abline") +
    xlab("RPKM 1 nucosol") + ylab("RPKM 2 nucosol") 
  ggsave(getEPlotPath("nuc-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
  
  ggplot(subset(df, biotype=="lnc" & (RPKM1.nuc < 500) & (RPKM2.nuc < 500)), 
         aes(x = RPKM1.nuc,y=RPKM2.nuc))  + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("nucleus replicates lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("RPKM 1 nucosol") + ylab("RPKM 2 nucosol") 
  ggsave(getEPlotPath("nuc-cellVbiotypeCompare-lncRNA-lncRNA-lt500.pdf.pdf"),height=24,width=7)
  
  
  ## RPKMratio -- cyt/nuc 
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)")
  ggsave(getEPlotPath("cytNucRatio-logDensity-biotype.pdf"),height=5,width=10)
  
  ggplot(df, aes(RPKMratio,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)") + xlim(0,20)
  ggsave(getEPlotPath("cytNucRatio-Density-biotype-lt20.pdf"),height=5,width=10)
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw()
  ggsave(getEPlotPath("cytNucRatio-logDensity-genetype.pdf"),height=5,width=10)
  
  ggplot(df, aes(RPKMratio,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) + xlim(0,20) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)") 
  ggsave(getEPlotPath("cytNucRatio-Density-genetype.pdf"),height=5,width=12)
  
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)")
  ggsave(getEPlotPath("cytNucRatio-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()
  ggsave(getEPlotPath("cytNucRatio-logDensity-genetypeVcell.pdf"),height=24,width=7)
 
  

  ### Cytosol Expression
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-biotype.pdf"),height=5,width=10)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-genetype.pdf"),height=2,width=10)
  
  ggplot(subset(df,RPKMsum.cyt < 20), aes(RPKMsum.cyt,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-Density-biotype-lt500.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.cyt < 20), aes(RPKMsum.cyt,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-Density-genetype-lt500.pdf"),height=5,width=12)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()+ 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-genetypeVcell.pdf"),height=24,width=7)
  
  ### Nuclear Expression
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-biotype.pdf"),height=7,width=12)
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-genetype.pdf"),height=7,width=15)
  
  ggplot(subset(df,RPKMsum.nuc < 20), aes(RPKMsum.nuc,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-Density-biotype-lt500.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.nuc < 20), aes(RPKMsum.nuc,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-Density-genetype-lt500.pdf"),height=5,width=12)
  
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()+ 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-genetypeVcell.pdf"),height=24,width=7)
  
}





cytNucRatiosPseudoCount <- function(){
  if(!file.exists(local.plotdir.encodeExpr)){
    dir.create(local.plotdir.encodeExpr,recursive=T)
  }
  df <- cytNucRatios()
  
  ### Get number of found genes/celltype
  df$Count <- 1
  dgc.table <- tbl_df(df)
  dgc.table.cell <- group_by(dgc.table, cell,gene_type)
  df.cell.counts <- as.data.frame(summarise(dgc.table.cell,
                                            num.types = n()))
  
  
  
  ggplot(df.cell.counts, aes(x=cell, y = num.types, fill=gene_type))+ geom_bar(stat="identity") +
    coord_flip() + theme_bw() + 
    ggtitle("genes w/ expr > 0 in (cyt,nuc)")
  ggsave(getEPlotPath("cellTypeCountByGeneType.pdf"),height=7,width=7) #1
  
  
  lncwFnc <- filter(df,biotype == "lnc") %.%
    group_by(gene_type,cell,funcLnc) %.%
    summarise(num.types = n()) 
  
  ggplot(as.data.frame(lncwFnc), aes(x=cell, y = num.types, fill=gene_type))+ geom_bar(stat="identity") +
    theme_bw() + facet_grid(funcLnc ~ . ,scale="free")+
    ggtitle("lncRNA w/ expr > 0 in (cyt,nuc)\n{0,1} == {func,unlabelled}") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(getEPlotPath("cellTypeCountByLncFunc.pdf"),height=7,width=7) #1
  
  
  
  
  funcLnc.biotype <- unique(filter(df, biotype == "lnc", funcLnc == 1) %.%
                              select(gene_id,gene_type)) %.%
    group_by(gene_type) %.%
    summarise(num.types = n())
  
  ggplot(funcLnc.biotype, aes(x=gene_type, y = num.types)) + geom_bar(stat="identity") +
    coord_flip() + theme_bw() +
    ggtitle("functional Lnc found by biotype")
  ggsave(getEPlotPath("funcLncCountByGeneType.pdf"),height=7,width=7) # 2
  
  ggplot(select(df,funcLnc==1), aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("functional Lnc found by biotype") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")
  ggsave(getEPlotPath("cytVsNuc-lncFunc.pdf"),height=7,width=7) # 3
  
  ggplot(subset(df,biotype == "lnc"), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc)))+
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("All Lnc: cyt vs. nuc\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")+
    scale_size_manual(values = c(3,4))
  ggsave(getEPlotPath("cytVsNuc-All-lnc.pdf"),height=7,width=7) # 4
  
  ggplot(subset(df,biotype == "lnc"& (RPKMsum.cyt < 500) & (RPKMsum.nuc < 500)), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc)))+
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    ggtitle("All Lnc: cyt vs. nuc\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum")+
    scale_size_manual(values = c(3,4))
  ggsave(getEPlotPath("cytVsNuc-All-lnc-lt500.pdf"),height=7,width=7) 
  
  
  
  ggplot(subset(df,biotype=="lnc"), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~gene_type, nrow=1,scale="free") +
    ggtitle("biotype=Lnc by gene_type\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,4))
  ggsave(getEPlotPath("cytVsNuc-lncBy.pdf"),height=5,width=10) # 5
  
  ggplot(arrange(subset(df,biotype=="lnc" & RPKMsum.nuc < 500 & RPKMsum.cyt < 500),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~gene_type, nrow=1,scale="free") +
    ggtitle("biotype=Lnc by gene_type\ny=x abline,RPKMsum... < 500") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3.5))
  ggsave(getEPlotPath("cytVsNuc-lncBy-lt500.pdf"),height=5,width=10) 
  
  ggplot(df,aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~biotype, nrow=1,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-biotypeCompare.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.nuc < 500 & RPKMsum.cyt < 500),aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point(alpha=I(0.5)) +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_wrap(~biotype, nrow=1,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-biotypeCompare-lt500.pdf"),height=5,width=10)
  
  
  ggplot(df,aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~biotype,scale="free") +
    ggtitle(" genes facet by biotype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare.pdf"),height=24,width=7)
  
  ## cyt/nuc sum comparison
  ggplot(subset(df, biotype=="lnc"), aes(x = RPKMsum.nuc,y=RPKMsum.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") 
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
  
  ggplot(arrange(subset(df, biotype=="lnc" & (RPKMsum.cyt < 500) & (RPKMsum.nuc < 500)),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3))
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA-lt500.pdf"),height=24,width=7)
  
  ggplot(arrange(subset(df, biotype=="lnc" & (RPKMsum.cyt < 100) & (RPKMsum.nuc < 100)),funcLnc), 
         aes(x = RPKMsum.nuc,y=RPKMsum.cyt,color=factor(funcLnc),size=factor(funcLnc))) + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cyt/nuc sum lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("nuclear RPKM sum") + ylab("cytosol RPKM sum") +
    scale_size_manual(values = c(2.4,3)) +xlim(0,100) + ylim(0,100)
  ggsave(getEPlotPath("cytVsNuc-cellVbiotypeCompare-lncRNA-lt100.pdf"),height=24,width=7)
  
  
  ## Correlation between replicates
  cyt.rsq <- as.data.frame(group_by(df,cell) %.% summarise(count = n()))
  cyt.rsq$count <- NULL
  cyt.rsq$gene_type <- "total"
  cyt.rsq$rsquaredCyt <- as.numeric(group_by(df,cell) %.% do(function(x)summary(lm(RPKM1.cyt ~ RPKM2.cyt, data=x))$r.squared))
  cyt.rsq$rsquaredNuc <- as.numeric(group_by(df,cell) %.% do(function(x)summary(lm(RPKM1.nuc ~ RPKM2.nuc, data=x))$r.squared))
  cyt.rsq.lnc <- as.data.frame(group_by(df,cell,gene_type) %.% summarise(count = n()))
  cyt.rsq.lnc$count <- NULL
  cyt.rsq.lnc$rsquaredCyt <- as.numeric(group_by(df,cell,gene_type) %.% do(function(x)summary(lm(RPKM1.cyt ~ RPKM2.cyt, data=x))$r.squared))
  cyt.rsq.lnc$rsquaredNuc <- as.numeric(group_by(df,cell,gene_type) %.% do(function(x)summary(lm(RPKM1.nuc ~ RPKM2.nuc, data=x))$r.squared))
  cytNuc.comb.rsq <- rbind(cyt.rsq.lnc,cyt.rsq)
  
  ggplot(melt(cytNuc.comb.rsq, id.var=c("cell","gene_type")), aes(x=cell,y=value,fill=value)) + geom_bar(stat="identity") +
    facet_grid(variable~gene_type,scale="free") + theme_bw() +
    xlab("R-squared for subset/total in cell type") + ylab("cell type") + 
    ggtitle("lm on RPKM1 & RPKM2 of cytosol\nlncRNA")  +coord_flip() +
    scale_fill_gradient2(low="red",mid="grey",high="green",midpoint=0.75)
  ggsave(getEPlotPath("cyt-Stats-allgenes.pdf"),height=7,width=7)
  
  
  
  ## cyt sum comparison
  
  
  ggplot(subset(df, biotype=="lnc"), 
         aes(x = RPKM1.cyt,y=RPKM2.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cytosols replicates \n lncRNA facet by genetype\ny=x abline") +
    xlab("RPKM 1 cytosol") + ylab("RPKM 2 cytosol") 
  ggsave(getEPlotPath("cyt-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
  
  ggplot(subset(df, biotype=="lnc" & (RPKM1.cyt < 500) & (RPKM2.cyt < 500)), 
         aes(x = RPKM1.cyt,y=RPKM2.cyt)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("cytosols replicates \n lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("RPKM 1 cytosol") + ylab("RPKM 2 cytosol") 
  ggsave(getEPlotPath("cyt-cellVbiotypeCompare-lncRNA-lncRNA-lt500.pdf.pdf"),height=24,width=7)
  
  
  ## nuc sum comparison
  ggplot(subset(df, biotype=="lnc"), 
         aes(x = RPKM1.nuc,y=RPKM2.nuc)) + 
    geom_point() +
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("nucleus replicates \n lncRNA facet by genetype\ny=x abline") +
    xlab("RPKM 1 nucosol") + ylab("RPKM 2 nucosol") 
  ggsave(getEPlotPath("nuc-cellVbiotypeCompare-lncRNA.pdf"),height=24,width=7)
  
  ggplot(subset(df, biotype=="lnc" & (RPKM1.nuc < 500) & (RPKM2.nuc < 500)), 
         aes(x = RPKM1.nuc,y=RPKM2.nuc))  + 
    geom_point()+
    geom_abline(slope=1, intercept=0) + theme_bw()+
    facet_grid(cell~gene_type,scale="free") +
    ggtitle("nucleus replicates lncRNA facet by genetype(0,500)\ny=x abline") +
    xlab("RPKM 1 nucosol") + ylab("RPKM 2 nucosol") 
  ggsave(getEPlotPath("nuc-cellVbiotypeCompare-lncRNA-lncRNA-lt500.pdf.pdf"),height=24,width=7)
  
  
  ## RPKMratio -- cyt/nuc 
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)")
  ggsave(getEPlotPath("cytNucRatio-logDensity-biotype.pdf"),height=5,width=10)
  
  ggplot(df, aes(RPKMratio,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)") + xlim(0,20)
  ggsave(getEPlotPath("cytNucRatio-Density-biotype-lt20.pdf"),height=5,width=10)
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw()
  ggsave(getEPlotPath("cytNucRatio-logDensity-genetype.pdf"),height=5,width=10)
  
  ggplot(df, aes(RPKMratio,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) + xlim(0,20) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)") 
  ggsave(getEPlotPath("cytNucRatio-Density-genetype.pdf"),height=5,width=12)
  
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("Ratio=(RPKMsum.cyt/RPKMsum.nuc)")
  ggsave(getEPlotPath("cytNucRatio-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMratio),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()
  ggsave(getEPlotPath("cytNucRatio-logDensity-genetypeVcell.pdf"),height=24,width=7)
  
  
  
  ### Cytosol Expression
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-biotype.pdf"),height=5,width=10)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-genetype.pdf"),height=2,width=10)
  
  ggplot(subset(df,RPKMsum.cyt < 20), aes(RPKMsum.cyt,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-Density-biotype-lt500.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.cyt < 20), aes(RPKMsum.cyt,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-Density-genetype-lt500.pdf"),height=5,width=12)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMsum.cyt),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()+ 
    ggtitle("SUM=(RPKM1.cyt + RPKM2.cyt)")
  ggsave(getEPlotPath("cytSum-logDensity-genetypeVcell.pdf"),height=24,width=7)
  
  ### Nuclear Expression
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-biotype.pdf"),height=7,width=12)
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-genetype.pdf"),height=7,width=15)
  
  ggplot(subset(df,RPKMsum.nuc < 20), aes(RPKMsum.nuc,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-Density-biotype-lt500.pdf"),height=5,width=10)
  
  ggplot(subset(df,RPKMsum.nuc < 20), aes(RPKMsum.nuc,fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(. ~ gene_type,scales="free") + theme_bw() + 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-Density-genetype-lt500.pdf"),height=5,width=12)
  
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ biotype,scales="free") + theme_bw() +
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-biotypeVcell.pdf"),height=24,width=7)
  
  ggplot(df, aes(log(RPKMsum.nuc),fill=factor(funcLnc))) + 
    geom_density(alpha=I(0.4)) +
    facet_grid(cell ~ gene_type,scales="free") + theme_bw()+ 
    ggtitle("SUM=(RPKM1.nuc + RPKM2.nuc)")
  ggsave(getEPlotPath("nucSum-logDensity-genetypeVcell.pdf"),height=24,width=7)
  
}










