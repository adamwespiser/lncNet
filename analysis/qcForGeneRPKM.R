







getReportsForGene <- function(){
  list(flux=readInTable(getFullPath("/data/fluxCapData-lpa-proc-REPORT.tab")),
       rpkmFromBam=readInTable(getFullPath("/data/rpkmFromBAMCapData-lpa-proc-REPORT.tab")),
       rsem=readInTable(getFullPath("/data/rsemCapData-lpa-proc-REPORT.tab")))
}

plotRepQC <- function(){
  rep = getReportsForGene()
  colNames <- c("cell", "localization", "replicate","genesFound", "experiment")
  flux=rep$flux[colNames]
  rsem=rep$rsem[colNames]
  rpkmFromBam=rep$rpkmFromBam[colNames]
  flux$method = "flux"
  rsem$method = "rsem"
  rpkmFromBam$method = "rpkmFromBam"
  comb <- rbind(flux,rsem,rpkmFromBam)
  comb$tag <- with(comb,paste(cell,experiment,sep="."))

  ggplot(comb, aes(x=tag, y=genesFound))+facet_grid(method~.) + geom_bar() +
    ggtitle("Number of genes with read information\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("count of genes found") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/foundGenes.pdf"),height=7,width=12)


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








