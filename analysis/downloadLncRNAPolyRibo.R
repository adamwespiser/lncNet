


metaURL <- "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB5049&result=read_run&fields="
metaLocal <- getFullPath("data/lncRNA_polyRibosome.tab")

dataDir <- "/project/umw_zhiping_weng/wespisea/lncPolyRibo/"

downloadMetaData <- function(){
  if(!file.exists(metaLocal)){
    download.file(url=metaURL,destfile=metaLocal)
  }
  meta.df <- read.csv(file=metaLocal,stringsAsFactors=FALSE,sep="\t")
  meta.df$localFile <- paste0(dataDir,meta.df$run_accession)
  meta.df$downloadCmd <- paste0("wget --continue ",meta.df$fastq_ftp," -O ", dataDir,meta.df$run_accession,".fastq.gz")
  
  cmd.vec <- meta.df$downloadCmd
  
    md <- paste0(cmd.vec,rep(c(" &"," "),length.out =length(cmd.vec)))
  write(md, file="~/sandbox/downloadLncPolyRibo.sh")
  scpFile(file.local="~/sandbox/downloadLncPolyRibo.sh", dir.remote="~/bin/")
  # chmod u+x ~/sandbox/downloadLncPolyRibo.sh;~/sandbox/downloadLncPolyRibo.sh
  
}