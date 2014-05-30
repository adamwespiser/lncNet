
filesTxtTab <<- "~/data/wgEncodeCshlLongRnaSeqFiles.tab"
rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"

genStarAlignCmdSpikeIN.paramFileUniq <- function(rd1,rd2,outfile){
  paramFile <- "/home/aw30w/log/params/parametersPolyRibo_AWmod.txt"
  paste0("STAR --runMode alignReads ", 
         " --readFilesIn ", rd1, " ", rd2, 
         " --outFileNamePrefix ", outfile,
         " --parametersFiles ", paramFile)
}


generateRPKMFromBamFromStar <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  # df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  #  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$starAln <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star.samAligned.out.sam"))
  df$overlapSize <- floor(as.numeric(gsub(x=gsub(x=df$readType,pattern="2x",replacement=""),pattern="D",replacement=""))/2)
  df.comb$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".uniq.star_sort.bam"))
  df.comb$bamBai <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam.bai"))
  
  
  cmd0 <- "grep @ test.star.samAligned.out.sam > test.uniq.star.samAligned.out.sam;;grep NH:i:1 test.star.samAligned.out.sam >> test.uniq.star.samAligned.out.sam"
  o0 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd0,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  cmd1 <- "samtools view -bS test.uniq.star.samAligned.out.sam -o test.uniq.star.bam;;samtools sort -m 171798691840 test.uniq.star.bam test.uniq.star_sort"
  o1 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  cmd2 <- "samtools index test.uniq.star_sort.bam"
  o2 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd2,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  
  cmd3 <- "java -jar -Xmx24g /home/aw30w/bin/bam2rpkm-0.06/bam2rpkm-0.06.jar -f /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.NIST14SpikeIn.gtf -i test.uniq.star_sort.bam --overlap xxOOxx -r exon -o test.uniq.transByExon.gtf"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  o3 <- as.character(unlist(sapply(seq_along(o3), function(i)gsub(x=o3[i],pattern="xxOOxx", replacement=df$overlapSize[i]))))
  
  
  outputTotal <- sapply(paste(o0,o1,o2,o3,sep=";;"), function(x)gsub(x=x,pattern="//",replacement="/"))
  # 183840
  
  
  write(outputTotal, file="~/sandbox/procUniqReads")
  scpFile(file.local="~/sandbox/procUniqReads", dir.remote="~/bin/")
  # cat ~/bin/procUniqReads | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 183840 -W 600 -Q short -t procUniqReads -i "{}"
  fileOut <- paste0(file.path(rnaseqdir,"starSpikeIn",df.comb[which(df.comb$read1.rnaExtract == "longPolyA"),"bare"]),".uniq.transByExon.gtf")
  sapply(fileOut, hpc.file.exists)                      
  
  
  fileOut <- paste0(file.path(rnaseqdir,"starSpikeIn",df.comb[which(df.comb$read1.rnaExtract == "longPolyA"),"bare"]),".uniq.star.samAligned.out.sam")
  sapply(fileOut, hpc.file.exists)  
  
}






generateRPKMFromBamFromSortedSam <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  # df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  #  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$starAln <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star.samAligned.out.sam"))
  df$overlapSize <- floor(as.numeric(gsub(x=gsub(x=df$readType,pattern="2x",replacement=""),pattern="D",replacement=""))/2)
  df.comb$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".uniq.star_sort.bam"))
  df.comb$bamBai <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam.bai"))
  
  
  cmd0 <- "head -n 50000 test.star_sort.sam | grep \\^@  > test.uniq.star_sort.sam;;grep NH:i:1 test.star_sort.sam >> test.uniq.star_sort.sam"
  o0 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd0,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  cmd1 <- "samtools view -bS  test.uniq.star_sort.sam -o test.uniq.star_sort.bam"
  o1 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  cmd2 <- "samtools index test.uniq.star_sort.bam"
  o2 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd2,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  
  
  cmd3 <- "java -jar -Xmx24g /home/aw30w/bin/bam2rpkm-0.06/bam2rpkm-0.06.jar -f /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.NIST14SpikeIn.gtf -i test.uniq.star_sort.bam --overlap xxOOxx -r exon -o test.uniq.star_sort.transByExon.gtf"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  o3 <- as.character(unlist(sapply(seq_along(o3), function(i)gsub(x=o3[i],pattern="xxOOxx", replacement=df$overlapSize[i]))))
  
  
  outputTotal <- sapply(paste(o0,o1,o2,o3,sep=";;"), function(x)gsub(x=x,pattern="//",replacement="/"))
  # 183840
  
  
  write(outputTotal, file="~/sandbox/procStarSort")
  scpFile(file.local="~/sandbox/procStarSort", dir.remote="~/bin/")
  # cat ~/bin/procStarSort | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 8 -m 18340 -W 600 -Q short -t starSort4 -i "{}"
  
  fileOut <- paste0(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),".uniq.star_sort.transByExon.gtf")
  endFileFound <- sapply(fileOut, hpc.file.exists)   
  
  fileOut <- paste0(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),".star_sort.sam")
  sapply(fileOut, hpc.file.exists)  
  
}
generateRPKMFromBamFromSortedSam <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  # df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  #  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  df$overlapSize <- floor(as.numeric(gsub(x=gsub(x=df$readType,pattern="2x",replacement=""),pattern="D",replacement=""))/2)
  
  
  
  o <- genStarAlignCmdSpikeIN.paramFileUniq(rd1=gsub(file.path(rnaseqdir,df.comb$read1.filename),pattern="\\.gz",replacement=""  ),
                                        rd2=gsub(file.path(rnaseqdir,df.comb$read2.filename),pattern="\\.gz",replacement=""  ),
                                        outfile=file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star.sam")))
  write(o, file="~/sandbox/runStar.sh")
  scpFile(file.local="~/sandbox/runStar.sh", dir.remote="~/bin/")
  # cat ~/bin/runStar.sh | xargs -I{}  perl ~/bin/runJob.pl -c 16 -m 3072 -W 600 -Q short -t "runStar" -i "{}"
  df.comb$starAln <- file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star.samAligned.out.sam"))
  
  
  cmd1 <- "samtools view -bS test.star.samAligned.out.sam -o test.star.bam;;samtools sort -m 171798691840 test.star.bam test.star_sort"
  o1 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn-uniq",filename)))))
  #write(o1, file="~/sandbox/starBamSamtools.sh")
  #scpFile(file.local="~/sandbox/starBamSamtools.sh", dir.remote="~/bin/")
  df.comb$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star_sort.bam"))
  cmd3 <- "samtools index test.star_sort.bam"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn-uniq",filename)))))
  write(o3, file="~/sandbox/stIdxStar")
  #scpFile(file.local="~/sandbox/stIdxStar", dir.remote="~/bin/")
  df.comb$bamBai <- file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star_sort.bam.bai"))
 
  cmd4 <- "java -jar -Xmx24g /home/aw30w/bin/bam2rpkm-0.06/bam2rpkm-0.06.jar -f /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.gtf -i test.star_sort.bam --overlap xxOOxx -r exon -o test.uniq.star_sort.transByExon.gtf"
  o4 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd4,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn-uniq",filename)))))
  o4 <- as.character(unlist(sapply(seq_along(o4), function(i)gsub(x=o4[i],pattern="xxOOxx", replacement=df$overlapSize[i]))))
  
  outputTotal <- as.character(unlist(sapply(paste(o,o1,o3,o4,sep=";;"), function(x)gsub(x=x,pattern="//",replacement="/"))))
  # 183840
  
  write(outputTotal, file="~/sandbox/starUniqr")
  scpFile(file.local="~/sandbox/starUniqr", dir.remote="~/bin/")
  # cat ~/bin/starUniqr | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 16 -m 6500 -W 600 -Q short -t starUniqr -i "{}"
  fileOut <-  file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star_sort.bam"))
  existFileout <- sapply(fileOut, hpc.file.exists)                      
  
  fileOut <-  file.path(rnaseqdir,"starSpikeIn-uniq",paste0(df.comb$bare,".star.samAligned.out.sam"))
  existFileout <- sapply(fileOut, hpc.file.exists)  
  outputTotal <- as.character(unlist(sapply(paste(o,o1,o3,sep=";;"), function(x)gsub(x=x,pattern="//",replacement="/"))))
  oNeed <- outputTotal[which(existsFileout == FALSE)]
  write(outputTotal, file="~/sandbox/su_miss")
  scpFile(file.local="~/sandbox/su_miss", dir.remote="~/bin/")
  # cat ~/bin/su_miss | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 16 -m 12500 -W 600 -Q short -t starUniqr -i "{}"
  
  
  
  cmd4 <- "java -jar -Xmx24g /home/aw30w/bin/bam2rpkm-0.06/bam2rpkm-0.06.jar -f /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.gtf -i test.star_sort.bam --overlap xxOOxx -r exon -o test.uniq.star_sort.transByExon.gtf"
  o4 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd4,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn-uniq",filename)))))
  o4 <- as.character(unlist(sapply(seq_along(o4), function(i)gsub(x=o4[i],pattern="xxOOxx", replacement=df$overlapSize[i]))))
  fileOuto4 <-  file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,"uniq.star_sort.transByExon.gtf"))
  existFileout <- sapply(fileOut, hpc.file.exists)  
  
  write(o4, file="~/sandbox/su_miss")
  scpFile(file.local="~/sandbox/su_miss", dir.remote="~/bin/")
  
  
  
  
  
  
  
  
  
  
}
