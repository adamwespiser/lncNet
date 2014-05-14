#! /usr/bin/Rscript --vanilla


readInGTExAllMeta <- function(){
  df <- read.csv(file = "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/GTExSraDB-metainfo.tab",stringsAsFactors=FALSE,sep="\t")
  datadir <- "/data/wespisea/gtex/fastq/"
  df$read1 <- paste0(datadir,df$run_accession,"_1.fasta.gz")
  df$read2 <- paste0(datadir,df$run_accession,"_2.fasta.gz")
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$downloadCmd <- paste0("/home/wespisea/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df
}


downloadFileMissing_url_getkey <- function(){
  df <- readInGTExAllMeta()
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$haveFiles <- file.exists(df$fastq1)
  df.need <- df[which(df$haveFiles == FALSE),]
  
  write(paste0(df.need$downloadCmd,rep(c(" &"," "),length=length(df.need$downloadCmd))), file="~/sandbox/downloadGTEx_need.sh") 
  # http://gap-upload.ncbi.nlm.nih.gov/E2004057-252A-4BC9-ADA3-2463E4AC9B98/SRR612551.sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47
  
  genWebsiteKey() #takes 30 seconds....
  key <- getWebsiteKey()
  print(paste("website key is ", key))
  o <- paste0("wget -O /data/wespisea/gtex/sra/",df.need$run_accession,".sra 'http://gap-upload.ncbi.nlm.nih.gov/",key,"/",
              df.need$run_accession,".sra?tic=ECA2A0FB-F5AC-43F1-BF95-7E38649C8A47' --continue",
              rep(c(" &"," "),length=length(df.need$run_accession)))
  
  len <- length(o)
  midpoint <- floor(len/2)
  write(o[1:midpoint], file = "~/sandbox/downloadGTEx2_wget1.sh")
  write(o[(midpoint+1):len], file = "~/sandbox/downloadGTEx2_wget2.sh")
  
  
  s0 <- "find /data/wespisea/gtex/sra/ -name \"*.sra\" -size 0 -delete"
  s1 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget1.sh\""
  s2 <- "screen -d -m sh -c \"~/sandbox/downloadGTEx2_wget2.sh\""
  cat(paste0(c(s0,s1,s2),collapse="\n"))
  c(s0,s1,s2)
  
}


runDownloadMoniter <- function(){
  
  while(TRUE){
    if (FALSE == SRAstatusGood()){
      cmds <- downloadFileMissing_url_getkey()
      system(cmds[1])
      Sys.sleep(10)
      system(cmds[2])
      system(cmds[3])
      Sys.sleep(10)
    } 
    Sys.sleep(3600)
  }
} 
