homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}


gtex.annot <<- getDataFile("GTExSampleAnnot.tab")
gtex.gene.expr <<- getDataFile("GTExGeneRPKM.gct")
gtex.gene.head <<- getDataFile("GTExGeneRPKM.gct")
k562.nuc.pap <<- getDataFile("wgEncodeCshlLongRnaSeqK562NucleusPapTranscriptGencV7.gtf")
k562.nuc.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqK562NucleusPapTranscriptGencV7.tab")
k562.cyt.pap <<- getDataFile("wgEncodeCshlLongRnaSeqK562CytosolPapTranscriptGencV7.gtf")
k562.cyt.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqK562CytosolPapTranscriptGencV7.tab")
IMR.cyt.pap <<-     getDataFile("wgEncodeCshlLongRnaSeqImr90CytosolPapTranscriptGencV10.gtf")
IMR.cyt.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqImr90CytosolPapTranscriptGencV10.tab")
IMR.nuc.pap <<-     getDataFile("wgEncodeCshlLongRnaSeqImr90NucleusPapTranscriptGencV10.gtf")
IMR.nuc.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqImr90NucleusPapTranscriptGencV10.tab")

convertENCODEGtfToTab <- function(gtfFile=k562.cyt.pap,tabFile=k562.cyt.pap.tab){
  system( paste("cat",gtfFile," | sed 's/[;\"]//g' |awk -F' ' '{print $2,$1,$4,$5,$10,$12,$6,$14,$16,$18}' > ",tabFile))
  
}

convertENCODEGtfsToTab <- function(){
  convertGtfToTab(k562.nuc.pap,k562.nuc.pap.tab)
  convertGtfToTab(k562.cyt.pap,k562.cyt.pap.tab)
  convertGtfToTab(IMR.cyt.pap,IMR.cyt.pap.tab)
  convertGtfToTab(IMR.nuc.pap,IMR.nuc.pap.tab)
}

readInENCODEGtfToTab <- function(file=IMR.nuc.pap.tab){
  ccol <- c("character","character", "numeric", "numeric","character", "character", rep("numeric",4))
  df <- read.csv(file=file, sep=" ", stringsAsFactors=FALSE,
                colClasses=ccol)
  colnames(df) <- c("version", "chr", "start", "stop", "gene_id", "transcript_id", "COMB", "RPKK1",
                    "RPKM2", "IDR")
}


getGTExAnnot <- function(file=gtex.annot){
  df <- read.csv(file,sep="\t",stringsAsFactors=FALSE)
  df$SAMPID <- as.character(sapply(df$SAMPID, function(x)gsub(x=x,pattern="[\\_\\-]",replacement="\\.")))
  df
}

getGTExAnnotExpr <- function(file=gtex.annot){
  df <- read.csv(file,sep="\t",stringsAsFactors=FALSE)
  df$SAMPID <- as.character(sapply(df$SAMPID, function(x)gsub(x=x,pattern="[\\_\\-]",replacement="\\.")))
  
  df <- df[grep("RNA", df$SMNABTCHT),]
  df <- df[which(df$SAMPID %in% getGTExGeneExprCols()),]
  df
}

getGTExGeneExprCols <- function(file=gtex.gene.expr){
  names(read.csv(file=gtex.gene.head, sep="\t",nrows=2,stringsAsFactors=FALSE))[c(-1,-2)]
}

getAnnotColNames <- function(file=gtex.annot){
  df <- read.csv(file=gtex.annot,sep="\t")
  colnames(df)
}

# cols can be either c("SMTS", "SMTSD")
getAnnotColTable <- function(file=gtex.annot,col="SMTS"){
  df <- getGTExAnnotExpr()
  table(df[[col]])
}

# example c("Nerve", "SMTS")
getGTExColsByAnnot <- function(query,colName="SMTS"){
  annot.df <- getGTExAnnot()

  if(!colName %in% colnames(annot.df)){
    stop("Colname not found in annotations")
  }
  all.cols <- names(table(annot.df[[colName]]))
  if(!(query %in% all.cols)){
    stop("colname not found")
  }
  else{
    annot.df$SAMPID[which(annot.df[[colName]] == query)]
  }
}

#expr.df <- getGeneExprForColMatch("Cells - Leukemia cell line (CML)", "SMTSD")
getGeneExprForColMatch <- function(query,colName="SMTS"){
  expr.cols <- getGTExGeneExprCols()
  colFound <- expr.cols %in% getGTExColsByAnnot(query,colName=colName)
  readCols <- c(NA,NA,ifelse(colFound == TRUE,NA,'NULL'))
  read.csv(file=gtex.gene.expr,stringsAsFactors=FALSE,colClasses=readCols,sep="\t")
}


convertToDf <- function(str,recSep,labelSep,valRm=""){
  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
 # vals <- sapply(vals, function(x)gsub(x,";",""))
  l <- as.list(vals)
  names(l) <- nms 
  as.data.frame(l)
}

shortenIdVec <- function(gene_id, sep = "\\."){
  as.character(sapply(gene_id, function(x)unlist(strsplit(x,split=sep))[1]))
}

convertGtfFile <- function(file=k562.cyt.pap,outfile=k562.cyt.pap.tab){
  
  d.ens <- read.csv(file=file,sep="\t",header=FALSE,
                colClasses=(c(rep('NULL',8), "character")),stringsAsFactors=FALSE)
  df <- do.call(rbind,lapply(d.ens$V9,function(x)convertToDf(x,"; ", " ",";")))
  
  
  df$iIDR <- as.numeric(sapply(levels(df$iIDR)[as.numeric(df$iIDR)], function(x)gsub(x,pattern=";",replacement="")))
  colClasses <- c("character", "character", "numeric", "numeric")
  for(i in seq_along(colClasses)){
    df[,i] <- eval(as.call(list(as.name(paste("as",colClasses[i], sep=".")), levels(df[,i])[as.numeric(df[,i])])))
  }
  df$gene_id_short <- shortenIdVec(df$gene_id)
  df$transcript_id_short <- shortenIdVec(df$transcript_id)
  
  d.first <- read.csv(file=f,sep="\t",header=FALSE,
                    colClasses=c(rep("character",3),rep("numeric",3), rep("character",1),"NULL","NULL"),stringsAsFactors=FALSE)
  colnames(d.first) <- c("chr","source","type","start","stop","COMB","strand")
  exportAsTable(file=outfile,cbind(d.first,df))
}
# file=k562.cyt.pap,outfile=k562.cyt.pap.tab




