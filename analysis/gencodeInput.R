gencodeV12.pc <<- getDataFile("gencodeV12.proteinCoding.gene.tab")
gencodeV12.lnc <<- getDataFile("gencodeV12.procTranLinc.gene.tab")

gencodeV10.lnc <<- getDataFile("gencode.v10.long_noncoding_RNAs.gtf")
gencodeV10.lnc.website <<- "ftp://ftp.sanger.ac.uk/pub/gencode/release_10/gencode.v10.long_noncoding_RNAs.gtf.gz"
gencodeV10.lnc.tab <<- getDataFile("gencode.v10.long_noncoding_RNAs.tab")
gencodeV10.pc.website <<- "ftp://ftp.sanger.ac.uk/pub/gencode/release_10/gencode.v10.pc_transcripts.fa.gz"
gencodeV10.pc <<- getDataFile("gencode.v10.pc_transcripts.fa")
gencodeV10.pc.tab <<- getDataFile("gencode.v10.pc_transcripts.tab")


gencodeV7.lnc <<- getDataFile("gencode.v7.long_noncoding_RNAs.gtf")
gencodeV7.lnc.website <<- "ftp://ftp.sanger.ac.uk/pub/gencode/release_7/gencode.v7.long_noncoding_RNAs.gtf.gz"
gencodeV7.lnc.tab <<- getDataFile("gencode.v7.long_noncoding_RNAs.tab")
gencodeV7.pc.website <<- "ftp://ftp.sanger.ac.uk/pub/gencode/release_7/gencode.v7.pc_transcripts.fa.gz"
gencodeV7.pc <<- getDataFile("gencode.v10.pc_transcripts.fa")
gencodeV7.pc.tab <<- getDataFile("gencode.v10.pc_transcripts.tab")


processGencodeV7.lnc <- function(file = gencodeV7.lnc, out = gencodeV7.lnc.tab){
  if(!file.exists(file)){
    download.file(url=gencodeV7.lnc.website,dest=paste(gencodeV7.lnc,".gz",sep=""))
    system(paste("gzip -d --force", paste(gencodeV7.lnc,".gz",sep="")))
  }
  
  p = pipe(paste("sed 's/[;\"]//g' ", file,  " |awk -F' ' '{ if ($3 == \"transcript\") print $10,$12,$14,$20}' |grep -E \"lincRNA|processed_transcript|antisense\""))
  df <- read.csv(file=p,sep=" ", stringsAsFactors=FALSE,blank.lines.skip=TRUE,header=FALSE) 
  exportAsTable(df=df,file=out)
}
getGencodeV7.lnc <- function(){
  if (!file.exists(gencodeV7.lnc.tab)){
    processGencodeV7.lnc()
  }
  df = read.csv(file=gencodeV7.lnc.tab,sep="\t", header=TRUE,stringsAsFactors=FALSE)
  colnames(df) <- c("gene_id", "transcript_id", "gene_type", "transcript_type")
  df
}
processGencodeV7.pc <- function(file = gencodeV10.pc, out = gencodeV10.pc.tab){
  if(!file.exists(file)){
    download.file(url=gencodeV10.pc.website,dest=paste(gencodeV10.pc,".gz",sep=""))
    system(paste("gzip -d --force", paste(gencodeV10.pc,".gz",sep="")))
  }
  
  p = pipe(paste("grep '>'", gencodeV10.pc ," | cut -d '|' -f 1,2  | sed 's/[>]//g' "))
  df <- read.csv(file=p,sep="|", stringsAsFactors=FALSE,blank.lines.skip=TRUE,header=FALSE) 
  colnames(df) <- c("transcript_id", "gene_id")
  exportAsTable(df=df,file=out)
}
getGencodeV7.pc <- function(){
  if (!file.exists(gencodeV7.pc.tab)){
    processGencodeV7.pc()
  }
  df = read.csv(file=gencodeV7.pc.tab,sep="\t",stringsAsFactors=FALSE)
  df
}



processGencodeV10.lnc <- function(file = gencodeV10.lnc, out = gencodeV10.lnc.tab){
  if(!file.exists(file)){
    download.file(url=gencodeV10.lnc.website,dest=paste(gencodeV10.lnc,".gz",sep=""))
    system(paste("gzip -d --force", paste(gencodeV10.lnc,".gz",sep="")))
  }
  
  p = pipe(paste("sed 's/[;\"]//g' ", file,  " |awk -F' ' '{ if ($3 == \"transcript\") print $10,$12,$14,$20}' |grep -E \"lincRNA|processed_transcript|antisense\""))
  df <- read.csv(file=p,sep=" ", stringsAsFactors=FALSE,blank.lines.skip=TRUE,header=FALSE) 
  exportAsTable(df=df,file=out)
}
getGencodeV10.lnc <- function(){
  if (!file.exists(gencodeV7.lnc.tab)){
    processGencodeV7.lnc()
  }
  df = read.csv(file=gencodeV10.lnc.tab,sep="\t", header=TRUE,stringsAsFactors=FALSE)
  colnames(df) <- c("gene_id", "transcript_id", "gene_type", "transcript_type")
  df
}

processGencodeV10.pc <- function(file = gencodeV10.pc, out = gencodeV10.pc.tab){
  if(!file.exists(file)){
    download.file(url=gencodeV10.pc.website,dest=paste(gencodeV10.pc,".gz",sep=""))
    system(paste("gzip -d --force", paste(gencodeV10.pc,".gz",sep="")))
  }
  
  p = pipe(paste("grep '>'", gencodeV10.pc ," | cut -d '|' -f 1,2  | sed 's/[>]//g' "))
  df <- read.csv(file=p,sep="|", stringsAsFactors=FALSE,blank.lines.skip=TRUE,header=FALSE) 
  colnames(df) <- c("transcript_id", "gene_id")
  exportAsTable(df=df,file=out)
}
getGencodeV10.pc <- function(){
  if (!file.exists(gencodeV10.pc.tab)){
    processGencodeV10.pc()
  }
  df = read.csv(file=gencodeV10.pc.tab,sep="\t", stringsAsFactors=FALSE)
  df
}



getLncGenesV12 <- function(file = gencodeV12.lnc){
  df <- read.csv(file=file, sep=" ", header=FALSE,colClass=c("character",'NULL',"character",'NULL'))
  colnames(df) <- c("gene_id","biotype")
  df$gene_id_short <- shortenIdVec(df$gene_id)
}

getPcGenesV12 <- function(file = gencodeV12.pc){
  df <- read.csv(file=file, sep=" ", header=FALSE,colClass=c("character",'NULL',"character",'NULL'))
  colnames(df) <- c("gene_id","biotype")
  df$gene_id_short <- shortenIdVec(df$gene_id)
  df
}

getGencodeAnnot <- function(biotype,gencodeVersion){
  if(!(biotype %in% c("lnc", "pc")) && !(gencodeVersion %in% c("v7", "v10")))
    stop("biotype must be either <lnc> or <pc> & gencodeVersion must be either <v7> or <v10>")
  df <- data.frame()
  if(identical(tolower(gencodeVersion), "v7")){
    if(identical(tolower(biotype), "lnc")){
      df <- getGencodeV7.lnc()
    } else if(identical(tolower(biotype), "pc")){
      df <- getGencodeV7.pc()
    }
  }
  else if (identical(tolower(gencodeVersion), "v10")){
    if(identical(tolower(biotype), "lnc")){
      df <- getGencodeV10.lnc()
    } else if(identical(tolower(biotype), "pc")){
      df <- getGencodeV10.pc()
    }
  } else if (identical(tolower(gencodeVersion), "v12")){
    if(identical(tolower(biotype), "lnc")){
      df <- getLncGenesV12()
    } else if(identical(tolower(biotype), "pc")){
      df <- getPcGenesV12()
    }
  }
  df
}
