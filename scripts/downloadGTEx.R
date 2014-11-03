#! /usr/bin/Rscript --vanilla
## USAGE:
#  run this program to ensure the GTEx download is continuously going
# screen -d -m sh -c "./downloadGTEx.R"

readInGTExAllMeta <- function(){
  df <- read.csv(file = "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/GTExSraDB-metainfo.tab",stringsAsFactors=FALSE,sep="\t")
  datadir <- "/data/wespisea/gtex/fastq/"
  df$read1 <- paste0(datadir,df$run_accession,"_1.fasta.gz")
  df$read2 <- paste0(datadir,df$run_accession,"_2.fasta.gz")
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$downloadCmd <- paste0("/home/wespisea/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df
}

genWebsiteKey <- function(SRA="SRR613771"){
  system(paste("cd /data/wespisea/gtex/sraDB; ~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/test-sra",SRA,"2>&1 > ~/sandbox/sratoolTest &"))
  Sys.sleep(30)
}
getWebsiteKey <- function(){
  lines <- readLines(con="~/sandbox/sratoolTest")
  lines <- lines[grep(x=lines, "Remote http")]
  lines <- lines[grep(x=lines, "http://gap-upload.ncbi.nlm.nih.gov/")]
  as.character(unlist(strsplit(x=as.character(unlist(strsplit(x=lines[1],split ="gap-upload.ncbi.nlm.nih.gov/")))[-1],split="/")[1]))[1]
}

downloadFileMissing_url_getkey <- function(){
  df <- readInGTExAllMeta()
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
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

SRAstatusGood <- function(sraDir="/data/wespisea/gtex/sra/"){
 
  sraDir <- "/data/wespisea/gtex/sra/"
  # sraDir <- "/home/wespisea/sandbox/testSRA"
  p <- pipe(paste("ls",sraDir))
  files <- file.path(sraDir,readLines(p))
  close(p)
  sizeDistro <- sapply(files, function(x)file.info(x)$size)
  if (sum(sizeDistro == 0) > 10){
    return(FALSE)
  }
  TRUE
}

runDownloadMoniter <- function(){
  print("entering download moniter")  
  i <- 1
  while(TRUE){
    print(paste("entering round:",i))
    status <- SRAstatusGood()
    print(paste("status is Good(t/f):",status))
    if (FALSE == SRAstatusGood()){
      print("download has failed, restarting")
      cmds <- downloadFileMissing_url_getkey()
      print("got new address, and commands")
      system(cmds[1])
      print("cleared size 0 files")
      Sys.sleep(10)
      print("entering new download commands to screen sessions")
      system(cmds[2])
      system(cmds[3])
      Sys.sleep(10)
    }
    Sys.sleep(600)
  }
} 


#runDownloadMoniter()

runConvertToFastq <- function(script="~/sandbox/downloadGTEx_sraToFastQ.sh"){
  #downloadFileMissing()
  df <- readInGTExAllMeta()
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip /data/wespisea/gtex/sra/",df$run_accession,".sra")
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
  df$haveFiles <- file.exists(df$fastq1)
  df.need <- df[which(df$haveFiles == FALSE),]
  
  
  o.convert <- paste0(df.need$cddownloadCmd,rep(c(" &"," "),length=length(df.need$cddownloadCmd)))
  
  len <- length(o.convert)
  midpoint <- floor(len/2)
  write(o.convert[1:midpoint], file = "~/sandbox/downloadGTEx_sraToFastQ_1.sh")
  write(o.convert[(midpoint+1):len], file = "~/sandbox/downloadGTEx_sraToFastQ_2.sh")
  
  convert.cmd1 <- paste0("screen -d -m sh -c \"~/sandbox/downloadGTEx_sraToFastQ_1.sh\"")
  convert.cmd2 <- paste0("screen -d -m sh -c \"~/sandbox/downloadGTEx_sraToFastQ_1.sh\"")
  system(convert.cmd1)
  system(convert.cmd2)
  
  
  # ../sra/SRR612347.sra 
}


readInGTExAllMeta <- function(){
  df <- read.csv(file = "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/GTExSraDB-metainfo.tab",stringsAsFactors=FALSE,sep="\t")
  datadir <- "/data/wespisea/gtex/fastq/"
  df$read1 <- paste0(datadir,df$run_accession,"_1.fasta.gz")
  df$read2 <- paste0(datadir,df$run_accession,"_2.fasta.gz")
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$downloadCmd <- paste0("/home/wespisea/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip ",df$run_accession)
  df$gtexId <- as.character(unlist(sapply(df$sample_attribute, function(x)convertToList(x)[["submitted_sample_id"]])))
  
  df
}
getGTExAnnot <- function(file= "/home/wespisea/data/GTExSampleAnnot.tab"){
  df <- read.csv(file,sep="\t",stringsAsFactors=FALSE)
  df$SAMPID <- as.character(sapply(df$SAMPID, function(x)gsub(x=x,pattern="[\\_\\-]",replacement="\\.")))
  df
}

getGTExWithFullAnnot <- function(){
  d <- getGTExAnnot()
  m <- readInGTExAllMeta()
  d$SAMPIDdash = gsub(d$SAMPID, pattern="\\.",replacement="-")
  m$gtexId= sapply(m$gtexId,toupper)
  d$SAMPIDdash= sapply(d$SAMPIDdash,toupper)
  comb <- merge(d,m, by.x="SAMPIDdash",by.y="gtexId")
  
}

getGTExWithMeta <- function(){
  annot <- getGTExWithFullAnnot()
 
  df <- readInGTExAllMeta()
  df$cddownloadCmd <- paste0("cd /data/wespisea/gtex/sraDB;~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -O /data/wespisea/gtex/fastq/ --split-files -gzip /data/wespisea/gtex/sra/",
                             df$run_accession,".sra")
  gtexDir <- "/data/wespisea/gtex/fastq/"
  df$fastq1 <- paste0(gtexDir, df$run_accession, "_1.fastq.gz")
  df$fastq2 <- paste0(gtexDir, df$run_accession, "_2.fastq.gz")
  df$sraFile <- paste0("/data/wespisea/gtex/sra/",df$run_accession,".sra")
  
  df$haveMetaData <- ifelse(df$run_accession %in% annot$run_accession, TRUE, FALSE)
  df$fastq1Exists <- sapply(df$fastq1, file.exists)
  df$fastq2Exists <- sapply(df$fastq2, file.exists)
  df  
}



runConvertToFastqMeta <- function(script="~/sandbox/downloadGTEx_sraToFastQ.sh"){
  #downloadFileMissing()
  df <- getGTExWithMeta()
  df.need <- df[-which(df$fastq1Exists  & df$fastq2Exists),]
  df.need.meta <- df.need[which(df.need$haveMetaData),]
  df.need.noMeta <- df.need[-which(df.need$haveMetaData),]
  
  #split need meta
  o.convert <- paste0(df.need.meta$cddownloadCmd,
                      rep(c(" &"," "),length=length(df.need.meta$cddownloadCmd)))
  
  len <- length(o.convert)
  midpoint <- floor(len/2)
  
  write(o.convert[1:midpoint], file = "~/sandbox/downloadGTEx_sraToFastQ_1.sh")
  write(o.convert[(midpoint+1):len], file = "~/sandbox/downloadGTEx_sraToFastQ_2.sh")
  write(df.need.noMeta$cddownloadCmd, file = "~/sandbox/downloadGTEx_sraToFastQ_3.sh")
  
  convert.cmd1 <- paste0("screen -d -m sh -c \"~/sandbox/downloadGTEx_sraToFastQ_1.sh\"")
  convert.cmd2 <- paste0("screen -d -m sh -c \"~/sandbox/downloadGTEx_sraToFastQ_2.sh\"")
  convert.cmd3 <- paste0("screen -d -m sh -c \"~/sandbox/downloadGTEx_sraToFastQ_3.sh\"")
  
  system(convert.cmd1)
  system(convert.cmd2)
  system(convert.cmd3)
  
  
  # ../sra/SRR612347.sra 
}

#/home/wespisea/bin/FastQC/fastqc -q -o ~/sandbox SRR1092397_1.fastq.gz

fastqcReport <- function(){
  df <- getGTExWithMeta()
  reportdir <- "/data/wespisea/gtex/fastq-FastQC-report/"
  df$fastqc1File <- paste0(reportdir, df$run_accession, "_1")
  df$fastqc2File <- paste0(reportdir, df$run_accession, "_2")
  
  df$cmd1fastqc <- paste0("mkdir ",df$fastqc1File, ";/home/wespisea/bin/FastQC/fastqc -q -o ",df$fastqc1File,"/ ", df$fastq1)
  df$cmd2fastqc <- paste0("mkdir ",df$fastqc2File,";/home/wespisea/bin/FastQC/fastqc -q -o ",df$fastqc2File,"/ ", df$fastq2)
  #/data/wespisea/gtex/fastq-FastQC-report/SRR615008_1/SRR615008_1_fastqc.html
  df$cmd1fastqcDest <- paste0("/data/wespisea/gtex/fastq-FastQC-report/",
                              df$run_accession,"_1/",df$run_accession,"_1_fastqc.html")
  df$cmd2fastqcDest <- paste0("/data/wespisea/gtex/fastq-FastQC-report/",
                              df$run_accession,"_2/",df$run_accession,"_2_fastqc.html")
  df$fastqc1Exists <- sapply(df$cmd1fastqcDest,file.exists)
  df$fastqc2Exists <- sapply(df$cmd2fastqcDest,file.exists)
  
  
  df.need <- df[which(df$fastq1Exists == FALSE & df$fastq2Exists == FALSE),]
  df.need.meta <- df.need[which(df.need$haveMetaData),]
  df.need.noMeta <- df.need[which(df.need$haveMetaData == FALSE),]

  df.have <- df[which(df$fastq1Exists  & df$fastq2Exists),]
  df.have.meta <- df.have[which(df.have$haveMetaData),]
  
  fastq1cmds <- df.have.meta$cmd1fastqc[which(df.have.meta$fastqc1Exists == FALSE)]
  fastq2cmds <- df.have.meta$cmd2fastqc[which(df.have.meta$fastqc2Exists == FALSE)]
  
  
  
  idx <- order(c(seq_along(fastq1cmds), seq_along(fastq2cmds)))
  cmds <- c(fastq1cmds,fastq2cmds)[idx]
  write(cmds, file = "~/sandbox/makeGTExFastqcReport.sh")
  fastqcCmd <- paste0("screen -d -m sh -c \"~/sandbox/makeGTExFastqcReport.sh\"")
  cat(fastqcCmd)
  
}

starGenerateGenome_ZLAB <- function(){
  rnaseqdir <- "/data/wespisea/STAR-genome/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDir/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  # STAR --runMode genomeGenerate --genomeDir genomepath --genomeFastaFiles  genomepath/genome.fa  
  # --sjdbGTFfile genomepath/genes.gtf --sjdbOverhang 75 --runThreadN 8
  generateGenome <- paste0("/home/wespisea/bin/STAR_2.3.0e.Linux_x86_64/STAR --runMode genomeGenerate --genomeDir ", starGenomeDir, " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf, " --sjdbOverhand 75 --runThreadN 6") 
  createGenome <- paste0("screen -d -m sh -c \"",generateGenome,"\"")
  cat(createGenome)
}

