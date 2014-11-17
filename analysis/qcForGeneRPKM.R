







getReportsForGene <- function(){
  list(flux=readInTable(getFullPath("/data/fluxCapData-lpa-proc-REPORT.tab")),
       rpkmFromBam=readInTable(getFullPath("/data/rpkmFromBAMCapData-lpa-proc-REPORT.tab")),
       rsem=readInTable(getFullPath("/data/rsemCapData-lpa-proc-REPORT.tab")),
       eXpress=readInTable(getFullPath("/data/eXpressCapData-lpa-proc-REPORT.tab")),
       rfgTopTransEx=readInTable(getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-RRPM-REPORT.tab")),
       rfgExUniqReads=readInTable(getFullPath("/data/rpkmFromBam-ExonCounting-TopTransCellType-UNIQ-RRPM-REPORT.tab")))}

plotRepQC <- function(){
  rep = getReportsForGene()
  colNames <- c("cell", "localization", "replicate","genesFound", "experiment")
  flux=rep$flux[colNames]
  rsem=rep$rsem[colNames]
  eXpress=rep$eXpress[colNames]
  
  rpkmFromBam=rep$rpkmFromBam[colNames]
  flux$method = "flux"
  rsem$method = "rsem"
  eXpress$method="eXpress"
  rpkmFromBam$method = "rpkmFromBam"
  rfgExUniqReads$metod = "rfgExUniqReads"
  comb <- rbind(flux,rsem,rpkmFromBam,rfgExUniqReads)
  comb$tag <- with(comb,paste(cell,experiment,sep="."))

  ggplot(comb, aes(x=tag, y=genesFound))+facet_grid(method~.) + geom_bar() +
    ggtitle("Number of genes with read information\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("count of genes found") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/foundGenes.pdf"),height=7,width=12)


}


plotRepQC_meanRPKM <- function(){
  rep = getReportsForGene()
  colNames <- c("cell", "localization","experiment", "replicate","genesExpressed", "meanRPKM","sumRPKM")
  
  flux=rep$flux[c("cell", "localization","experiment", "replicate","genesExpressed", "meanRPKMforGene")]
  flux$sumRPKM = with(flux,meanRPKMforGene * genesExpressed)
  rpkmFromBam=rep$rpkmFromBam[c("cell", "localization","experiment", "replicate","genesExpressed", "meanRPKM","sumRPKM")]
  rfgTopTransEx=rep$rfgTopTransEx[c("cell", "localization","experiment", "replicate","genesExpressed", "meanRPKM","sumRPKM")]
  rfgExUniqReads=rep$rfgExUniqReads[c("cell", "localization","experiment", "replicate","genesExpressed", "meanRPKM","sumRPKM")]
  
  rsem=rep$rsem[c("cell", "localization","experiment", "replicate","genesExpressed", "meanFPKM","sumFPKM")]
  eXpress=rep$eXpress[c("cell", "localization","experiment", "replicate","genesExpressed", "meanFPKM","sumFPKM")]
  
  
  colnames(flux) <- colNames
  colnames(rsem) <- colNames
  
  
  flux$method = "flux"
  rsem$method = "rsem"
  eXpress$method = "eXpress"
  
  rpkmFromBam$method = "rpkmFromBam"
  rfgTopTransEx$method = "rfgTopTransEx"
  rfgExUniqReads$method = "rfgExUniqReads"
  
  comb <- rbind(flux,rsem,rpkmFromBam,rfgTopTransEx,rfgExUniqReads)
  comb$tag <- with(comb,paste(cell,experiment,sep="."))
  
  ggplot(comb, aes(x=tag, y=genesExpressed))+facet_grid(method~.,scale="free_y") + geom_bar() +
    ggtitle("Number of genes with read information\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("count of genes expressed") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/genesExpressed.pdf"),height=9,width=12)
  
  ggplot(comb, aes(x=tag, y=meanRPKM))+facet_grid(method~.,scale="free_y") + geom_bar() +
    ggtitle("meanRPKM\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("meanRPKM") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/meanRPKM.pdf"),height=9,width=12)
  
  ggplot(comb, aes(x=tag, y=sumRPKM))+facet_grid(method~.,scale="free_y") + geom_bar() +
    ggtitle("sumRPKM\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("sumRPKM") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/sumRPKM.pdf"),height=9,width=12)
  
  
  
}

plotReadLengths <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df$readLength <- as.numeric(gsub(sapply(strsplit(df$readType, "x"), function(x)x[2]),pattern="D",replacement=""))
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")

#/project/umw_zhiping_weng/wespisea/rna-seq//starSpikeIn/
df.comb$remote <- file.path(rnaseqdir,"starSpikeIn/",paste0(df.comb$bare,".trans.gtf"))

# if(writeCopyScript){
#   o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, " /home/wespisea/data/rpkmFromBam/",paste0(df.comb$bare,".trans.gtf"))
#   write(o1,file="~/sandbox/rpkmFromBamFetch")
# }
df.comb$rpkmFromBamFile <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".trans.gtf")
df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"rpkmFromBamFile", "bare","read1.readLength")]
df.comb$rfbGene <- paste0("/home/wespisea/data/rpkmFromBam/",df.comb$bare,".genes.gtf")
colnames(df.comb) <- c("localization", "cell", "rnaExtract", "replicate", 
                       "rpkmFromBamFile", "bare", "readLength","rfbGene")
ggplot(subset(df.comb,rnaExtract=="longPolyA"), aes(x=replicate,y=readLength,fill=factor(readLength))) + 
  geom_bar(stat="identity") + 
  facet_grid(cell ~ localization)+
  ylim(0,101) + theme_bw() + 
  ggtitle("readLength of cell/localization/replicate\nRNA-seq data from CSHL")
ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/rnaSeq-readLength.pdf"),height=12,width=5)

}








