#source("dataInput.R")

getStarSpikeIn <- function(){
  df <- read.csv(file ="~/sandbox/countFound_exonMerge.stat",stringsAsFactors=FALSE,header=FALSE,sep=" ")
  colnames(df) <- c("count", "file")
 df$name <- as.character(unlist(sapply(df$file, function(x)unlist(strsplit(x,split="/")[1])[7])))
 df$base <- as.character(unlist(sapply(df$name, function(x)unlist(strsplit(x,split="\\.")[1])[1])))
  df$gene_region <- as.character(unlist(sapply(df$name, function(x)unlist(strsplit(x,split="\\.")[1])[2])))
  
  dd <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  dd <- subset(dd, type == "bam" &  (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  dd$name <- as.character(unlist(sapply(dd$filename, function(x)unlist(strsplit(x,split="\\.")[1])[1])))
  dd$name <- as.character(unlist(sapply(dd$name, function(x)gsub(x, pattern="Aln", replacement="Fastq"))))
  
   merge.df <- merge(df,dd, by.x="base", by.y="name")
   m.df <- merge.df[c("cell","localization","rnaExtract","count","gene_region")]
   m.df$expr <- paste0(m.df$localization,".",m.df$rnaExtract)
   
   sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  
   ggplot(sum.df[which(sum.df$gene_region != "gene"),], aes(x=expr,y=count,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap("cell") + coord_flip()
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/totalReads-combined.pdf"))
  
  sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  sum.df <- sum.df[which(sum.df$gene_region != "gene"),]
  s.df <- ddply(sum.df,.(expr,cell), transform, sumCount = sum(count) )
  
  ggplot(s.df, aes(x=expr,y=count/sumCount,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap("cell") + coord_flip()
  
  m.df <- merge.df[c("cell","localization","rnaExtract","count","gene_region","replicate")]
  m.df$expr <- paste0(m.df$localization,".",m.df$rnaExtract,".",m.df$replicate)
  m.df <- m.df[which(m.df$gene_region != "gene"),]
  s.df <- ddply(m.df,.(expr,cell), transform, sumCount = sum(count) )
  

  ggplot(s.df[order(s.df$gene_region),],  aes(x=expr,y=count/sumCount,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap(~cell) + coord_flip() + ggtitle("STAR mapped to gencodeV19 introns/exons")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/STAR-fracReadsMapped-byReplicate.pdf"))
  
  ggplot(s.df[order(s.df$gene_region),],  aes(x=expr,y=count,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap(~cell) + coord_flip() + ggtitle("STAR mapped to gencodeV19 introns/exons")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/STAR-totalReadsMapped-byReplicate.pdf"))
  
  sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  
}

getSegMappedRegions <- function(){
  df <- read.csv(file ="~/sandbox/segCountFound.space",stringsAsFactors=FALSE,header=FALSE,sep=" ")
  colnames(df) <- c("count", "file")
  df$name <- df$file
  #df$name <- as.character(unlist(sapply(df$file, function(x)unlist(strsplit(x,split="/")[1])[7])))
  df$base <- as.character(unlist(sapply(df$file, function(x)unlist(strsplit(x,split="\\.")[1])[1])))
  df$gene_region <- as.character(unlist(sapply(df$name, function(x)unlist(strsplit(x,split="\\.")[1])[2])))
  
  dd <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  dd <- subset(dd, type == "bam" &  (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  dd$name <- as.character(unlist(sapply(dd$filename, function(x)unlist(strsplit(x,split="\\.")[1])[1])))
  dd$name <- as.character(unlist(sapply(dd$name, function(x)gsub(x, pattern="Aln", replacement="Fastq"))))
  
  merge.df <- merge(df,dd, by.x="base", by.y="name")
  m.df <- merge.df[c("cell","localization","rnaExtract","count","gene_region")]
  m.df$expr <- paste0(m.df$localization,".",m.df$rnaExtract)
  
  sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  
  ggplot(sum.df[which(sum.df$gene_region != "gene"),], aes(x=expr,y=count,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap("cell") + coord_flip()
  ggsave(getFullPath("plots/rnaExpr/mappedReads/segemehl/totalReads-combined.pdf"))
  
  sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  sum.df <- sum.df[which(sum.df$gene_region != "gene"),]
  s.df <- ddply(sum.df,.(expr,cell), transform, sumCount = sum(count) )
  
   
  m.df <- merge.df[c("cell","localization","rnaExtract","count","gene_region","replicate")]
  m.df$expr <- paste0(m.df$localization,".",m.df$rnaExtract,".",m.df$replicate)
  m.df <- m.df[which(m.df$gene_region != "gene"),]
  s.df <- ddply(m.df,.(expr,cell), transform, sumCount = sum(count) )
  
  
  ggplot(s.df[order(s.df$gene_region),],  aes(x=expr,y=count/sumCount,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap(~cell) + coord_flip() + ggtitle("segemehl mapped to gencodeV19 introns/exons")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/segemehl/segemehl-fracReadsMapped-byReplicate.pdf"))
  
  ggplot(s.df[order(s.df$gene_region),],  aes(x=expr,y=count,fill=gene_region))+geom_bar(stat="identity") +
    facet_wrap(~cell) + coord_flip() + ggtitle("segemehl mapped to gencodeV19 introns/exons")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/segemehl/segemehl-totalReadsMapped-byReplicate.pdf"))
  
  
  sum.df <- ddply(m.df, .(expr,cell,gene_region),summarise, count=sum(count))
  
}


