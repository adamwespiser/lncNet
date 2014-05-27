


metaURL <- "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB5049&result=read_run&fields="
metaLocal <- getFullPath("data/lncRNA_polyRibosome.tab")

dataDir <<- "/project/umw_zhiping_weng/wespisea/lncPolyRibo/"

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

genStarAlignCmdSpikeIN.paramFile.polyR <- function(rd1,outfile){
  paramFile <- "/home/aw30w/log/params/parametersPolyRibo_AWmod.txt"
  paste0("STAR --runMode alignReads ", 
         " --readFilesIn ", rd1,
         " --outFileNamePrefix ", outfile,
         " --parametersFiles ", paramFile)
}
#  /home/aw30w/bin/bwa-0.5.9/bwa aln -c -l 25 -k 2 -n 10 -a bwtsw -p /project/umw_zhiping_weng/wespisea/lncPolyRibo//bwa-index/GRCh37.p13.genome /project/umw_zhiping_weng/wespisea/lncPolyRibo/ERR380658.fastq.gz
# /home/aw30w/bin/bwa-0.5.9/bwa aln -c -l 25 -k 2 -n 10 /project/umw_zhiping_weng/wespisea/lncPolyRibo/bwa-index/GRCh37.p13.genome /project/umw_zhiping_weng/wespisea/lncPolyRibo/ERR380651.fastq.gz 


genBwaAlignCmdSpikeIN.paramFile.polyR <- function(rd1,outfile){
  paramFile <- "/project/umw_zhiping_weng/wespisea/lncPolyRibo//bwa-index/GRCh37.p13.genome"
  paste0("/home/aw30w/bin/bwa-0.5.9/bwa aln -c -l 25 -k 2 -n 10 ", 
         paramFile,
         " ", rd1,
         " > ", outfile)
}

generateRPKMFromBamFromStar_lncPoly <- function(){
 
  rnaseqdir <- dataDir
  
  df <- read.csv(file=metaLocal, stringsAsFactors=FALSE, sep="\t")
  # df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  #  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  meta.df$localFile <- paste0(dataDir,meta.df$run_accession)
  meta.df$read1 <- paste0(dataDir,meta.df$run_accession,".fastq.gz")
  
  meta.df$downloadCmd <- paste0("wget --continue ",meta.df$fastq_ftp," -O ", meta.df$read1)
  
  read1 <- meta.df$localFile  
  
  meta.df$bare <- meta.df$localFile
  meta.df$starAln <- paste0(file.path(rnaseqdir,"starSpikeIn",meta.df$run_accession),"uniq.star.samAligned.out.sam")
  meta.df$overlapSize <- 20
  meta.df$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn",paste0(meta.df$run_accession,".uniq.star_sort.bam"))
  meta.df$bamBai <- file.path(rnaseqdir,"starSpikeIn",paste0(meta.df$run_accession,".star_sort.bam.bai"))
  
  rnaSeqDirUltimate <- "bam059"
  cat("mkdir ",file.path(rnaseqdir,rnaSeqDirUltimate))
  
  
  o <- genBwaAlignCmdSpikeIN.paramFile.polyR(rd1=meta.df$read1,
                                        outfile=file.path(rnaseqdir,rnaSeqDirUltimate,paste0(meta.df$run_accession,".bam")))
  
  

 # cmd0 <- "grep @ test.star.samAligned.out.sam > test.uniq.star.samAligned.out.sam;;grep NH:i:1 test.star.samAligned.out.sam >> test.uniq.star.samAligned.out.sam"
#  o0 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd0,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  cmd1 <- "samtools sort -m 85798691840 test.bam test.bwa_sort"
  o1 <- as.character(unlist(sapply(meta.df$run_accession, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,rnaSeqDirUltimate,filename)))))
  
  cmd2 <- "samtools index test.bwa_sort.bam"
  o2 <- as.character(unlist(sapply(meta.df$run_accession, function(filename)gsub(x=cmd2,pattern="test", replacement=file.path(rnaseqdir,rnaSeqDirUltimate,filename)))))
  
  cmd3 <- "java -jar -Xmx24g /home/aw30w/bin/bam2rpkm-0.06/bam2rpkm-0.06.jar -f /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.NIST14SpikeIn.gtf -i test.bwa_sort.bam --overlap xxOOxx -r exon -o test.transByExon.gtf"
  o3 <- as.character(unlist(sapply(meta.df$run_accession, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,rnaSeqDirUltimate,filename)))))
  o3 <- as.character(unlist(sapply(seq_along(o3), function(i)gsub(x=o3[i],pattern="xxOOxx", replacement="20"))))
  
  outputTotal <- as.character(unlist(sapply(paste(o,o1,o2,o3,sep=";;"), function(x)gsub(x=x,pattern="//",replacement="/"))))
  # 183840
  
  write(outputTotal, file="~/sandbox/lncPolyStar")
  scpFile(file.local="~/sandbox/lncPolyStar", dir.remote="~/bin/")
  # cat ~/bin/lncPolyStar | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 16 -m 6500 -W 600 -Q short -t mapPolyLnc -i "{}"
  fileOut <- paste0(rnaseqdir,rnaSeqDirUltimate,"/",meta.df$run_accession,".bam")
  sapply(fileOut, hpc.file.exists)                      
  
  fileOut <- paste0(rnaseqdir,rnaSeqDirUltimate,"/",meta.df$run_accession,".transByExon.gtf")
  sapply(fileOut, hpc.file.exists)                      
  
  
}







