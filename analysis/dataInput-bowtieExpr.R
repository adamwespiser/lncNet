homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }

source(getFullPath("analysis/gencodeInput.R"))
source(getFullPath("analysis/clusterExecute.R"))

filesTxtTab <<- "~/data/wgEncodeCshlLongRnaSeqFiles.tab"
local.datadir <<- "/home/wespisea/data/"
rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"
 
ifNFileExists <- function(f,cmd){
  paste0("[[ ! -s ",f," ]] && ",cmd)
}

ifFileExists <- function(f,cmd){
  paste0("[[ -s ",f," ]] && ",cmd)
}

ifNBothFilesExist <- function(f1,f2,cmd){
  paste0("([[ ! -s ",f1," ]] || ","[[ ! -s ",f2," ]]) && ( ",cmd," )")
}

ifFileExistsDelete <- function(checkFile,delFile){
  paste0("( [[  -s ",checkFile," ]] && ","rm ( " ,delFile," ) )")
}

getRnaSeqFiles <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & rnaExtract == "longPolyA")
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  # rsem-calculate-expression --paired-end wgEncodeCshlLongRnaSeqSknshraCellLongnonpolyaFastqRd1Rep1.fastq.gz wgEncodeCshlLongRnaSeqSknshraCellPapFastqRd2Rep1.fastq.gz /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref-spikeIn/ -p 8 --ci-memory 8G --samtools-sort-mem 8GB 
  df.comb$read1 <- file.path(rnaseqdir,df.comb$read1.filename)
  df.comb$read2 <- file.path(rnaseqdir,df.comb$read2.filename)
  df.comb$read1.fa <- gsub(x=df.comb$read1,pattern=".gz",replacement="")
  df.comb$read2.fa <- gsub(x=df.comb$read2,pattern=".gz",replacement="")
  
  df.comb$rsemOutput <- file.path(rnaseqdir, "star-transcriptome","RSEM",df.comb$bare)
  df.comb$rsemOutputE2e <- file.path(rnaseqdir, "star-transcriptome","RSEM-e2e",df.comb$bare)
  
  df.comb$rpkmFromBamOutput <- file.path(rnaseqdir, "star-transcriptome","RSEM","rpkmFromBam",df.comb$bare)
  df.comb$starOutput <- file.path(rnaseqdir, "star-transcriptome",paste0(df.comb$bare,".star.sam"))
  df.comb$starAln <- paste0(starOutput, "Aligned.out.sam") # AlignedToTranscriptome.out.bam
  df.comb$AlignedToTrans <- paste0(starOutput, "AlignedToTranscriptome.out.bam") 
  
  df.comb$starOutputRSEM <- file.path(rnaseqdir, "star-transcriptome",paste0(df.comb$bare,".rsem.star.sam"))
  df.comb$starAlnRSEM <- paste0(df.comb$starOutputRSEM, "Aligned.out.sam") # AlignedToTranscriptome.out.bam
  df.comb$AlignedToTransRSEM <- paste0(df.comb$starOutputRSEM, "AlignedToTranscriptome.out.bam")  #Log.final.out
  df.comb$AlignedToTransRSEMFinalOut <- paste0(df.comb$starOutputRSEM, "Log.final.out")
  df.comb$rsemOutE2e <- file.path(rnaseqdir, "star-transcriptome","RSEM-e2e",df.comb$bare)
  df.comb$rsemOutFr <- file.path(rnaseqdir, "star-transcriptome","RSEM-fr",df.comb$bare)
  
  
  df.comb$starAlnBam <- paste0(starOutput, ".bam")
  df.comb$starSortBamShort <- paste0(starOutput, ".sort")
  df.comb$starSortBam <- paste0(starOutput, ".sort.bam")
  df.comb$expressOut <- file.path(rnaseqdir, "star-transcriptome","eXpress",df.comb$bare)
  df.comb$rsemOut <- file.path(rnaseqdir, "star-transcriptome","RSEM",df.comb$bare)
  
  df.comb$starExpressResults <-  file.path(rnaseqdir, "star-transcriptome","eXpress",df.comb$bare,"results.xprs")
  df.comb$starRsemResults <-  file.path(rnaseqdir, "star-transcriptome","RSEM-e2e",paste0(df.comb$bare,".genes.results"))
  
  df.comb$bowtieBam <- file.path(rnaseqdir, "bowtie",paste0(df.comb$bare,".transcript.bam"))
  df.comb$bowtieBamLong <- file.path(rnaseqdir, "bowtie","long",paste0(df.comb$bare,".transcript.bam"))
  
  
  df.comb$bowtieExpressOut <- file.path(rnaseqdir, "bowtie","eXpress",df.comb$bare)
  df.comb$bowtieRsemOut <- file.path(rnaseqdir, "bowtie","RSEM",df.comb$bare)
  
  df.comb$bowtieExpressResults <-  file.path(rnaseqdir, "bowtie","eXpress",df.comb$bare,"results.xprs")
  df.comb$bowtieRsemResults <-  file.path(rnaseqdir, "bowtie","RSEM",paste0(df.comb$bare,".genes.results"))
  
  df.comb$rpkmFromBamOutput <- file.path(rnaseqdir, "star-transcriptome","RSEM","rpkmFromBam",df.comb$bare)
  
  df.comb$StarRsemStats <- paste0(df.comb$AlignedToTransRSEM,".stats")
  df.comb$StarExpStats <- paste0(df.comb$AlignedToTrans,".stats")
  df.comb$bowtieLongStats <- paste0(df.comb$bowtieBamLong,".stats")
  
  
  
  df.comb
}
  
genStarAlignCmdTranscriptome <- function(rd1,rd2,outfile){
  #STAR --genomeDir /path/to/genome/ --readFilesIn Read1.gz Read2.gz --outSAMattributes NH   HI    
  #--outFilterMultimapNmax 20   --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04  
  #--alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8
  #--alignSJDBoverhangMin 1 --quantMode TranscriptomeSAM --runThreadN 12 &
  
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDirTrans-test/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  
  paste0("/home/aw30w/bin/STAR_2.3.1z4/STAR --runMode alignReads --genomeLoad LoadAndRemove --runThreadN 16", 
         " --readFilesIn ", rd1, " ", rd2, 
         " --genomeDir ", starGenomeDir,
         " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf,
         " --outFileNamePrefix ", outfile,
         " --outSAMattributes NH HI --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 ",
         " --alignIntronMin 2  --alignIntronMax 1 --alignMatesGapMax 1000000 --alignSJoverhangMin 8",
         " --alignSJDBoverhangMin 1 --quantMode TranscriptomeSAM ")
}

genStarAlignCmdTranscriptomeRSEM <- function(rd1,rd2,outfile){
  #STAR --genomeDir /path/to/genome/ --readFilesIn Read1.gz Read2.gz --outSAMattributes NH   HI    
  #--outFilterMultimapNmax 20   --outFilterMismatchNmax 999   --outFilterMismatchNoverLmax 0.04  
  #--alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8
  #--alignSJDBoverhangMin 1 --quantMode TranscriptomeSAM --runThreadN 12 &
  
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDirTrans-test/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  
  paste0("/home/aw30w/bin/STAR_2.3.1z4/STAR --runMode alignReads --genomeLoad LoadAndRemove --runThreadN 16", 
         " --readFilesIn ", rd1, " ", rd2, 
         " --genomeDir ", starGenomeDir,
         " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf,
         " --outFileNamePrefix ", outfile,
         " --outSAMattributes NH HI --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 ",
         " --alignIntronMin 2  --alignIntronMax 1 --alignMatesGapMax 1000000 --alignSJoverhangMin 8",
         " --alignSJDBoverhangMin 1 --quantMode TranscriptomeSAM ",
         " --alignEndsType EndToEnd ")
}

mapStarExpress <- function(){
  df.comb <- getRnaSeqFiles()
  o1 <- genStarAlignCmdTranscriptome(df.comb$read1.fa,df.comb$read2.fa,paste0(starOutput,".star.sam"))
  
  cmd4 <- paste0("/home/aw30w/bin/express-1.5.1-linux_x86_64/express  --no-update-check  ",
                 "--output-dir=", df.comb$expressOut,
                 " / /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19.transcripts.fa ",
                 df.comb$AlignedToTrans)
  o1.test <- ifNFileExists(df.comb$starExpressResults,o1)
  
  write(o1.test, file="~/sandbox/starTrans")
  scpFile(file.local="~/sandbox/starTrans", dir.remote="~/bin/")
  
  write(cmd4, file="~/sandbox/starTransExpress")
  scpFile(file.local="~/sandbox/starTransExpress", dir.remote="~/bin/")
  
  l1 <- 'cat ~/bin/starTrans | xargs -I{} perl /home/aw30w/bin/runJob.pl  -c 16 -m 3072 -W 240 -Q short -t "star1" -i "{}" | /home/aw30w/bin/getJobId > tmp111'
  l2 <- 'paste -d "#" tmp111  ~/bin/starTransExpress | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 5 -m 24576 -W 240 -Q short -t "starTransExpress" -i "{}" | /home/aw30w/bin/getJobId > tmp222'
  cat(l1,"\n",l2,"\n")

}

mapStarRSEM <- function(){
  df.comb <- getRnaSeqFiles()
  
  
  o1 <- genStarAlignCmdTranscriptomeRSEM(df.comb$read1.fa,df.comb$read2.fa,df.comb$starOutputRSEM)
  o1.test <- ifNFileExists(df.comb$AlignedToTransRSEMFinalOut,o1)
  o1.test1 <- ifNFileExists(df.comb$AlignedToTransRSEMFinalOut,paste0("echo ",o1))
  

  cmd2 <- paste0("rsem-calculate-expression -p 12 --paired-end --strand-specific --no-bam-output --bam  --num-threads 8 ",
                 df.comb$AlignedToTransRSEM, " /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19 ",
                 df.comb$rsemOutFr)
  cmd2.test1 <- ifNFileExists(df.comb$AlignedToTransRSEMFinalOut,paste0(" echo ",cmd2))
  
  
  
  cleanUp <- ifFileExistsDelete(checkFile =df.comb$AlignedToTransRSEMFinalOut,  df.comb$AlignedToTransRSEM)
  cleanUp.test1 <- ifFileExists(df.comb$AlignedToTransRSEMFinalOut,paste0(" rm  ",df.comb$AlignedToTransRSEM))
  
  write(o1.test1, file="~/sandbox/starRsemMapReads")
  scpFile(file.local="~/sandbox/starRsemMapReads", dir.remote="~/bin/")
  
  write(cmd2.test1, file="~/sandbox/starRsemQuant")
  scpFile(file.local="~/sandbox/starRsemQuant", dir.remote="~/bin/")
  
  write(cleanUp.test1, file="~/sandbox/starRsemClean")
  scpFile(file.local="~/sandbox/starRsemClean", dir.remote="~/bin/")
  
  
  l1 <- 'sh /home/aw30w/bin/starRsemMapReads | xargs -I{} perl /home/aw30w/bin/runJob.pl   -c 16 -m 4072 -W 1000 -Q long -t "srm_3" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/starRsem1'
  l2 <- 'paste -d "#" /home/aw30w/log/starRsem1  <(sh /home/aw30w//bin/starRsemQuant) | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 8 -m 10576 -W 240 -Q short -t "srq_3" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/starRsem2'
  l3 <- 'paste -d "#" /home/aw30w/log/starRsem2  <(sh /home/aw30w/bin/starRsemClean) | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 1 -m 1057 -W 240 -Q short -t "starRsemClean" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/starRsem3'
  
  cat(l1,"\n",l2)
  #  ifNFileExists(df.comb$AlignedToTransRSEMFinalOut,paste("echo, \"found\" "))

  #cat(l1,"\n",l2,"\n",l3,"\n")
  
}
  

genBowtieForStarComp <- function(rd1,rd2,outfile){
  bowtieRef <- "/project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19"
 #/share/pkg/bowtie/1.0.0/bowtie -q -n 2 -e 99999999 -l 25 -I 1 -X 1000 -p 16 -a -m 20 -S $rnaseqdata/rsem-ref/hg19gencodeV19 -1 read1.fastq -2 read2.fastq | /share/pkg/RSEM/1.2.11/sam/samtools view -S -b -o $rnaseqdata/rsem-hg19-gencodeV19-stranded/basefile.temp/basefile.bam
  
 paste0("/share/pkg/bowtie/1.0.0/bowtie  -n 2 -e 99999999 -l 25 -I 1 -X 1000 -p 16 -a -m 20 ",
          " -S ", bowtieRef, 
         " -1 ", rd1, " -2 ", rd2, 
        " | /share/pkg/RSEM/1.2.11/sam/samtools view -S -b -o ", outfile, " - ")
}
  
genBowtieForStarComp2 <- function(rd1,rd2,outfile,params,
                                  bowtieRef="/project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19"){
  #/share/pkg/bowtie/1.0.0/bowtie -q --e 99999999 -l 25 -I 1 -X 1000 -p 16 -a -m 20 -S $rnaseqdata/rsem-ref/hg19gencodeV19 -1 read1.fastq -2 read2.fastq | /share/pkg/RSEM/1.2.11/sam/samtools view -S -b -o $rnaseqdata/rsem-hg19-gencodeV19-stranded/basefile.temp/basefile.bam
  
  paste0("/share/pkg/bowtie/1.0.0/bowtie",params, " ",
         bowtieRef, 
         " -1 ", rd1, " -2 ", rd2, 
         " | /share/pkg/RSEM/1.2.11/sam/samtools view -S -b -o ", outfile, " - ")
}

mapBowtieRSEMcompExpressTEST <- function(){
  
  genFiles <- "
  cd $rnaseqdata
  mkdir -p test/RSEM
  mkdir -p test/eXpress
   head -n 50000 wgEncodeCshlLongRnaSeqHelas3CytosolPapFastqRd1Rep1.fastq > helaCytRd1Rep1.test.fastq
   head -n 50000 wgEncodeCshlLongRnaSeqHelas3CytosolPapFastqRd2Rep1.fastq > helaCytRd2Rep1.test.fastq
  "
  
  testDir <- "/project/umw_zhiping_weng/wespisea/rna-seq/test/"
  testRd1 <- paste0(testDir,"helaCytRd1Rep1.test.fastq")
  testRd2 <- paste0(testDir,"helaCytRd2Rep1.test.fastq")
  
  testBowtieResult <-   paste0(testDir,"helaCytRd2Rep1.transcript.bam")
  testRsemOut <-   paste0(testDir,"RSEM")
  testeXpressOut <-   paste0(testDir,"eXpress")
  
  test.bq <- genBowtieForStarComp(testRd1,testRd2,testBowtieResult)
  test.brq <- paste0("rsem-calculate-expression  --paired-end --no-bam-output --estimate-rspd --bam ",
                 testBowtieResult, " /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19 ",
                 testRsemOut)
  test.beq <- paste0("/home/aw30w/bin/express-1.5.1-linux_x86_64/express  --no-update-check --fr-stranded ",
                 "--output-dir=", testeXpressOut,
                 " /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19.transcripts.fa ",
                 testBowtieResult)
  cat(paste0(collapse="\n",c(genFiles,test.bq,test.brq,test.beq)))
  
  testBowtieResult2 <-   paste0(testDir,"helaCytRd2Rep1_v2.transcript.bam")
  testRsemOut2 <-   paste0(testDir,"RSEM2")
  testeXpressOut2 <-   paste0(testDir,"eXpress2")
  
  genDirs2 <- paste0(collapse=" && ",paste0("mkdir -p ",c(testRsemOut2,testeXpressOut2)))
  
  #   -X 800 -p 16 -k 20 --offrate 1
  test.bq2 <- genBowtieForStarComp2(testRd1,testRd2,testBowtieResult2,params=" -S -X 800 -p 16 -k 20 --offrate 1 ")
  test.brq2 <- paste0("rsem-calculate-expression -p 12 --paired-end --no-bam-output --estimate-rspd --bam ",
                     testBowtieResult2, " /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19 ",
                     testRsemOut2)
  test.beq2 <- paste0("/home/aw30w/bin/express-1.5.1-linux_x86_64/express  --no-update-check --fr-stranded ",
                     "--output-dir=", testeXpressOut2,
                     "  /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19.transcripts.fa ",
                    #" /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.bedtools.getfasta.split.fa ",
                     testBowtieResult2)
  #/project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.bedtools.getfasta.fa
  
  
  t2 <- paste0(collapse="\n",c(genDirs2,test.bq2,test.beq2,test.brq2))
  
  write(t2,file = "~/sandbox/bowtieTest")
  scpFile(file.local="~/sandbox/bowtieTest", dir.remote="~/bin/")
  
}




mapBowtieRSEMcompExpress <- function(long=FALSE,repl=FALSE,runNumber=2,INDEX="ALL"){
  df.comb <- getRnaSeqFiles()
  bowtieResult<- df.comb$bowtieBam
  
  
  if(identical(INDEX,"ALL")){
    jobsIndex=seq_along(df.comb$bowtieBam)
    
  } else {
    jobsIndex = INDEX
    print("only using indexes of: ")
    jobsIndex
  }
  
  df.comb <- df.comb[jobsIndex,]
  
  fileAppend <- character(length=0)
  if (TRUE == long){
    fileAppend <- "long"
    bowtieResult <- df.comb$bowtieBamLong
  } 
  # try ?best "-n 
  o1 <- genBowtieForStarComp2(df.comb$read1.fa,df.comb$read2.fa,bowtieResult,params=" -S -X 800 -p 16 -k 20 --offrate 1 ")
  o1.test <- ifNBothFilesExist(df.comb$bowtieExpressResults,df.comb$bowtieRsemResults,o1)
 
  cmd2 <- paste0("rsem-calculate-expression -p 12 --paired-end --strand-specific --no-bam-output --estimate-rspd --bam ",
                 bowtieResult, " /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19 ",
                 df.comb$bowtieRsemOut)
  cmd2.test <- ifFileExists(df.comb$bowtieBam,ifNFileExists(df.comb$bowtieRsemResults, cmd2))
  
  cmd3 <- paste0("/home/aw30w/bin/express-1.5.1-linux_x86_64/express  --no-update-check --fr-stranded ",
                 "--output-dir=", df.comb$bowtieExpressOut,
                 "  /project/umw_zhiping_weng/wespisea/rna-seq/rsem-ref/hg19gencodeV19.transcripts.fa ",
                 bowtieResult)
  cmd3.test <- ifFileExists(bowtieResult,ifNFileExists(df.comb$bowtieExpressResults, cmd3))
  
  cleanUp <- ifFileExists(bowtieResult, 
                          ifFileExists(df.comb$bowtieExpressResults, 
                                       paste("rm ", df.comb$bowtieBam)))
  
  hpc.file.exists(df.comb$bowtieRsemResults)
  
  write(o1, file=paste0("~/sandbox/bowtieCompareMap",fileAppend))
  scpFile(file.local=paste0("~/sandbox/bowtieCompareMap",fileAppend), dir.remote="~/bin/")
  
  write(cmd2, file=paste0("~/sandbox/bowtieCompareRSEM",fileAppend))
  scpFile(file.local=paste0("~/sandbox/bowtieCompareRSEM",fileAppend), dir.remote="~/bin/")

  write(cmd3, file=paste0("~/sandbox/bowtieCompareeXpress",fileAppend))
  scpFile(file.local=paste0("~/sandbox/bowtieCompareeXpress",fileAppend), dir.remote="~/bin/")
  
  write(cleanUp, file=paste0("~/sandbox/bowtieCompareCLEAN",fileAppend))
  scpFile(file.local=paste0("~/sandbox/bowtieCompareCLEAN",fileAppend), dir.remote="~/bin/")
  
  
#   l1 <- 'cat ~/bin/bowtieCompareMap | xargs -I{} perl /home/aw30w/bin/runJob.pl   -c 16 -m 12072 -W 240 -Q short -t "bowtieCompareMap" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp1'
#   l2 <- 'paste -d "#" /home/aw30w/log/bowtieExp1  ~/bin/bowtieCompareRSEM | xargs -I{} perl /home/aw30w/bin/runJobDep.pl  -c 12 -m 12576 -W 240 -Q short -t "bowtieCompareRSEM" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp2'
#   l3 <- 'paste -d "#" /home/aw30w/log/bowtieExp1  ~/bin/bowtieCompareeXpress | xargs -I{} perl /home/aw30w/bin/runJobDep.pl -c 5 -m 24576 -W 240 -Q short -t "bowtieCompareeXpress" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp3'
#   l4 <- 'paste -d "#" /home/aw30w/log/bowtieExp2  ~/bin/bowtieCompareCLEAN | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 1 -m 1057 -W 240 -Q short -t "bowtieCompareCLEAN" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp4'
  
  
  
  if (TRUE == long){
  l1 <- 'cat /home/aw30w/bin/bowtieCompareMap | xargs -I{} perl /home/aw30w/bin/runJob.pl   -c 16 -m 12072 -W 1440 -Q long -t "bmXXXX" -i "{}" | /home/aw30w/bin/getJobId |tee /home/aw30w/log/bowtieExp1_XXXX'
  l2 <- 'paste -d "#" /home/aw30w/log/bowtieExp1_XXXX  /home/aw30w/bin/bowtieCompareRSEM | xargs -I{} perl /home/aw30w/bin/runJobDep.pl  -c 12 -m 12576 -W 1440 -Q long -t "brqXXXX" -i "{}" | /home/aw30w/bin/getJobId | tee /home/aw30w/log/bowtieExp2_XXXX'
  l3 <- 'paste -d "#" /home/aw30w/log/bowtieExp1_XXXX  /home/aw30w/bin/bowtieCompareeXpress | xargs -I{} perl /home/aw30w/bin/runJobDep.pl -c 5 -m 10576 -W 1440 -Q long -t "beqXXXX" -i "{}" | /home/aw30w/bin/getJobId | tee /home/aw30w/log/bowtieExp3_XXXX'
  l4 <- 'paste -d "#" /home/aw30w/log/bowtieExp2_XXXX  /home/aw30w/bin/bowtieCompareCLEAN | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 1 -m 1057 -W 240 -Q short -t "bcXXXX" -i "{}" | /home/aw30w/bin/getJobId | tee /home/aw30w/log/bowtieExp4_XXXX'
  t.all <- paste(gsub(pattern="bowtieCompareMap",x=l1,replace=paste0("bowtieCompareMap",fileAppend)),
                gsub(pattern="bowtieCompareRSEM",x=l2,replace=paste0("bowtieCompareRSEM",fileAppend)),
                gsub(pattern="bowtieCompareeXpress",x=l3,replace=paste0("bowtieCompareeXpress",fileAppend)),
                 sep="\n")
  all.temp <- gsub(x=t.all,pattern="bowtieExp",replacement="bowtieExp_long")
  all <- gsub(x=all.temp,pattern="XXXX",replacement=runNumber)
  

  
  write(all, file=paste0("~/sandbox/bowtieLong",fileAppend))
  scpFile(file.local=paste0("~/sandbox/bowtieLong",fileAppend), dir.remote="~/bin/")
  
  if(TRUE == repl){
    print("sh /home/aw30w/bin/bowtieLonglong")
    cat(all) 
  }
  
  # end long
  } else {
    l1 <- 'cat ~/bin/bowtieCompareMap | xargs -I{} perl /home/aw30w/bin/runJob.pl   -c 16 -m 12072 -W 240 -Q short -t "bowtieCompareMap" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp1'
    l2 <- 'paste -d "#" /home/aw30w/log/bowtieExp1  /home/aw30w//bin/bowtieCompareRSEM | xargs -I{} perl /home/aw30w/bin/runJobDep.pl  -c 12 -m 12576 -W 240 -Q short -t "bowtieCompareRSEM" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp2'
    l3 <- 'paste -d "#" /home/aw30w/log/bowtieExp1  /home/aw30w/bin/bowtieCompareeXpress | xargs -I{} perl /home/aw30w/bin/runJobDep.pl -c 5 -m 24576 -W 240 -Q short -t "bowtieCompareeXpress" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp3'
    l4 <- 'paste -d "#" /home/aw30w/log/bowtieExp2  /home/aw30w/bin/bowtieCompareCLEAN | xargs -I{} perl /home/aw30w/bin/runJobDep.pl   -c 1 -m 1057 -W 240 -Q short -t "bowtieCompareCLEAN" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp4'
    out <- paste0(l1,"\n",l2,"\n",l3,"\n")
    
    
    write(out, file=paste0("~/sandbox/bowtieRun",fileAppend))
    scpFile(file.local=paste0("~/sandbox/bowtieRun",fileAppend), dir.remote="~/bin/")
    if(TRUE == repl){
      cat(l1,"\n",l2,"\n",l3,"\n")
    }
  }
}

getBamStats <- function(){
  
  df.comb <- getRnaSeqFiles()
  
  # df.comb$AlignedToTransRSEMFinalOut
  # df.comb$AlignedToTrans
  # df.comb$bowtieBamLong
  
  
  
  
}



#  [[ -s ~/ls.out ]] && echo "next" >> ~/ls.out && echo "next"
#  [[ -s file1 ]] && rm file2 && rm file3

#  [[ ! -s ~/ls.out ]] && [[ ! -s ~/jobs ]] && echo "next"
# ([[ -s a ]] || [[ -s b ]]) && echo "TRUE"
  

#bowtie for RSEM/exp done
#STAR -> express 
# express from bowtie needed
# RSEM from bowtie


# /home/aw30w/work/software/bam-parse/bin/test

getCounts <- function(){
  df.comb <- getRnaSeqFiles()
  inputfiles <- c(
  df.comb$AlignedToTransRSEM,
  df.comb$bowtieBamLong)
  outputfiles <- paste0(infiles, ".cnt")
  cmd1 <- paste0("/home/aw30w/work/software/bam-parse/bin/test --inputfile ",inputfiles,
                 " --outputfile ",outputfiles)
  write(cmd1, file="~/sandbox/bamCount")
  scpFile(file.local="~/sandbox/bamCount", dir.remote="~/bin/")
  
  
  
  l1 <- 'cat /home/aw30w/bin/bamCount | xargs -I{} perl /home/aw30w/bin/runJob.pl   -c 2 -m 4072 -W  240 -Q short -t "bc1" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/starRsem1'
  cat(l1)
  
  o1 <- genStarAlignCmdTranscriptome(df.comb$read1.fa,df.comb$read2.fa,paste0(starOutput,".star.sam"))
  input <- df.comb$AlignedToTrans
  output <- paste0(input,".cnt")
  cmd2 <- paste0("/home/aw30w/work/software/bam-parse/bin/test --inputfile ",input,
                 " --outputfile ",output)
  cmd3 <- paste0("rm ",df.comb$AlignedToTrans)
  write(paste(o1,cmd2,cmd3,sep=" && "), file="~/sandbox/star-bamCount")
  scpFile(file.local="~/sandbox/star-bamCount", dir.remote="~/bin/")
  l2 <- 'cat ~/bin/star-bamCount | xargs -I{} perl /home/aw30w/bin/runJob.pl  -c 16 -m 3072 -W 240 -Q short -t "statCount" -i "{}" | /home/aw30w/bin/getJobId > tmp111'
  cat(l1,"\n",l2)
}

jobs <- function(){
  df.comb <- getRnaSeqFiles()
  df.comb$AlignedToTransRSEMFinalOut
  
  srq <- 'cat  ~/bin/starRsemQuant | xargs -I{} perl /home/aw30w/bin/runJobRSEM.pl   -c 8 -m 10576 -W 1000 -Q long -t "srq_1" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/starRsem2'
  srq.output <- paste0("ls /project/umw_zhiping_weng/wespisea/rna-seq//star-transcriptome/RSEM-fr/*genes.results | wc -l")
  
  
  beq <- 'cat ~/bin/bowtieCompareeXpresslong | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 5 -m 10576 -W 1000 -Q long -t "beq" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp3'
  beq.output <- paste0("ls /project/umw_zhiping_weng/wespisea/rna-seq//bowtie/eXpress/*/results.xprs | wc -l ")
 # df.comb$bowtieExpressResults
  brq <- 'cat  ~/bin/bowtieCompareRSEMlong | xargs -I{} perl /home/aw30w/bin/runJobRSEM.pl  -c 12 -m 5128 -W 1000 -Q long -t "brq" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp2'
  brq.output <- paste0("ls /project/umw_zhiping_weng/wespisea/rna-seq//bowtie/RSEM-fr/*genes.results | wc -l ")
  #df.comb$bowtieRsemResults
  
  cat(srq,"\n",beq, "\n",brq)
 cat(srq.output,"\n",beq.output, "\n",brq.output)
 
 
 #963@ghpcc06:eXpress$ cat  ~/bin/bowtieCompareeXpresslong | sed 's:/long::g' > ~/bin/bowtieCompareeXpresslong2
 #964@ghpcc06:eXpress$ cat  ~/bin/bowtieCompareRSEMlong | sed 's:/long::g' > ~/bin/bowtieCompareRSEMlong2
 
 #
 beq.pre <- "cat  ~/bin/bowtieCompareeXpresslong | sed 's:/long::g' > ~/bin/bowtieCompareeXpresslong2"
 beq2 <- 'cat ~/bin/bowtieCompareeXpresslong2 | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 5 -m 10576 -W 1000 -Q long -t "beq2" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp3'
 beq.output <- paste0("ls /project/umw_zhiping_weng/wespisea/rna-seq//bowtie/eXpress/*/results.xprs | wc -l ")
 # df.comb$bowtieExpressResults
 brq.pre <- "cat  ~/bin/bowtieCompareRSEMlong | sed 's:/long::g' > ~/bin/bowtieCompareRSEMlong2"
 brq2 <- 'cat  ~/bin/bowtieCompareRSEMlong2 | xargs -I{} perl /home/aw30w/bin/runJobRSEM.pl  -c 12 -m 5128 -W 1000 -Q long -t "brq2" -i "{}" | /home/aw30w/bin/getJobId > /home/aw30w/log/bowtieExp2'
 brq.output <- paste0("ls /project/umw_zhiping_weng/wespisea/rna-seq//bowtie/RSEM/*genes.results | wc -l ")
 #df.comb$bowtieRsemResults
 
 cat(beq2, "\n",brq2)
 cat(beq.output, "\n",brq.output)
 
  
  runIndex <- which(hpc.file.exists(df.comb$bowtieRsemResults) == FALSE)
  mapBowtieRSEMcompExpress(long=TRUE,repl=TRUE,runNumber=4,INDEX=runIndex)
 
 df.comb$StarRsemStats <- paste0(df.comb$AlignedToTransRSEM,".stats")
 df.comb$StarExpStats <- paste0(df.comb$AlignedToTrans,".stats")
 df.comb$bowtieLongStats <- paste0(df.comb$bowtieBamLong,".stats")

}


  