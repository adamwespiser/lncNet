

thisTheme <<- theme_bw()+
  theme(text = element_text(size=9)) + 
  theme(panel.grid.major.y = element_line(colour = "grey"))+
  theme(panel.grid.minor.y = element_line(colour = "grey")) +
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank()) +
  theme(legend.position = "top")


plotCellTypeReadCountsCytNuc <- function(f=getFullPath("data/multicovCountsGraceFeatures.tab")){
  df <- read.csv(file=f,stringsAsFactors=FALSE,header=TRUE,sep="\t")
  df$counts <- as.numeric(df$counts)
  # colnames(df) <- c("cell","localization","rnaExtract","replicate","multicovGrace","cmd","cmdOut","counts")
  df$symbolShort <- paste0(ifelse(df$localization == "cytosol", "C","N"),".",ifelse(df$rnaExtract == "longNonPolyA","A-","A+"))
  df$symbol     <- paste0(  df$symbolShort,".R",df$replicate)
  symbolOrder <- c("N.A-","N.A+", "C.A+", "C.A-")
  df$order1 <- sapply(df$symbolShort, function(x)which(symbolOrder == x))
  df$order2 <- df$replicate
  df$cellTypeSample <- paste0(df$cell,df$symbol)
  
  celltypes = unique(df$cell)
  library(dplyr)
  
  df.symbol <- as.data.frame(group_by(df, cellTypeSample) %.%
                             mutate(total = sum(counts)) %.%
                             arrange(desc(order1),order2))
   
   plot.colors <- scale_fill_manual(values=c("cds"= rgb(0.3765,0.5843,0.7882),"utr3"= rgb(0.8,0.4,0.3725),"utr5"= rgb(0.6667, 0.7686,0.4235),
                                              "intron"= rgb(0.4824,0.098,0.4745), "noncoding"=rgb(0.3490,0.7294,0.8157)))
  
   for ( cell in celltypes){
     df.local <- df.symbol[which(df.symbol$cell == cell),]
     title.local <- paste0(cell," Mapped Reads\n","STAR aligned to hg19 & gencode.v19")
     file.base <- paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),cell)
     df.local$region <- factor(df.local$region, levels = rev(c("cds", "utr3", "utr5", "intron", "noncoding")))
     ggplot(df.local, aes(x=symbol,y=counts/total,fill=region)) + geom_bar(stat="bin") + 
       ggtitle(title.local) + plot.colors + thisTheme+ xlab("localization.pulldown.rep") + ylab("Fraction Uniquely Mapped Reads")+
       scale_x_discrete(limits=unique(as.data.frame(arrange(df.local,order1))[c("symbol")])$symbol) +
       theme(text = element_text(size=12))
     ggsave(file=paste0(file.base,"-frac-reads.pdf"),height=7,width=7)
     
     ggplot(df.local, aes(x=symbol,y=counts,fill=region)) + geom_bar(stat="bin") + 
       ggtitle(title.local) + plot.colors + thisTheme + xlab("localization.pulldown.rep") + ylab("Uniquely Mapped Reads")+
       scale_x_discrete(limits=unique(as.data.frame(arrange(df.local,order1))[c("symbol")])$symbol) +
       theme(text = element_text(size=12))
     ggsave(file=paste0(file.base,"-count-reads.pdf"),height=7,width=7)
     
   }
  
  df.symbol$region <- factor(df.symbol$region, levels = rev(c("cds", "utr3", "utr5", "intron", "noncoding")))
  df.symbol$symbol <- as.character(unlist(sapply(df.symbol$symbol, function(x)gsub(x=gsub(x,pattern="3",replace="1"),pattern="4",replacement="2"))))
  ggplot(df.symbol, aes(x=symbol,y=counts/total,fill=region)) + geom_bar(stat="bin") + 
    ggtitle("All celltypes: Mapped Reads") + plot.colors + thisTheme + xlab("localization.pulldown.rep") + ylab("Fraction Uniquely Mapped Reads")+
    scale_x_discrete(limits=unique(as.data.frame(arrange(df.symbol,order1))[c("symbol")])$symbol)  +
    facet_wrap(~cell,ncol=2,scale="free_x")
  ggsave(file=paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"allCells-frac-reads.pdf"),height=15,width=8)
  
  ggplot(df.symbol, aes(x=symbol,y=counts,fill=region)) + geom_bar(stat="bin",width=0.75) + 
    ggtitle("All celltypes: Mapped Reads") + plot.colors + thisTheme+ xlab("localization.pulldown.rep") + ylab("Uniquely Mapped Reads")+
    scale_x_discrete(limits=unique(as.data.frame(arrange(df.symbol,order1))[c("symbol")])$symbol)  +
    facet_wrap(~cell,ncol=2,scale="free_x") 

  ggsave(file=paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"allCells-total-reads.pdf"),height=15,width=8)
  
  
  
}