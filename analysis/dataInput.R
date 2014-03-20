homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}

local.datadir <<- "/home/wespisea/data/"

cshl.remove.rnaseq <<- c("wgEncodeCshlLongRnaSeqNhekCellPapTranscriptGencV7Rep5.gtf.gz","wgEncodeCshlLongRnaSeqNhekCellPamTranscriptGencV7Rep5.gtf.gz")

rnaexpr.comb.v7 <<-  getDataFile("cshl-rnaSeq-v7-allcombined.space")
rnaexpr.comb.v10 <<-  getDataFile("cshl-rnaSeq-v10-allcombined.space")

cshl.rnaseq.dir <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"

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
IMR.comb.expr   <<- getDataFile("ImrCytNucTopCombRnaExprTransByGene.tab")

fibroFile <<- getDataFile("transformedFibroblastsGTEx.tab")
gencode.pc <<- getDataFile("gencodeV12.proteinCoding.gene.tab")
gencode.lnc <<- getDataFile("gencodeV12.procTranLinc.gene.tab")

getLncGenesV12 <- function(file = gencode.lnc){
  df <- read.csv(file=file, sep=" ", header=FALSE,colClass=c("character",'NULL',"character",'NULL'))
  colnames(df) <- c("gene_id","biotype")
  df$gene_id_short <- shortenIdVec(df$gene_id)
  df
}

getPcGenesV12 <- function(file = gencode.pc){
  df <- read.csv(file=file, sep=" ", header=FALSE,colClass=c("character",'NULL',"character",'NULL'))
  colnames(df) <- c("gene_id","biotype")
  df$gene_id_short <- shortenIdVec(df$gene_id)
  df
}


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
  df <- read.csv(file=file, sep=" ", stringsAsFactors=FALSE, header=FALSE,
                 colClasses=ccol)
  colnames(df) <- c("version", "chr", "start", "stop", "gene_id", "transcript_id", "COMB", "RPKK1",
                    "RPKM2", "IDR")
  df
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
# getGTExColsByAnnot(query="Cells - Transformed fibroblasts",colName="SMTSD")
getGTExColsByAnnot <- function(query,colName="SMTS"){
  annot.df <- getGTExAnnotExpr()
  
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
# getGeneExprForColMatch(query="Cells - Transformed fibroblasts",colName="SMTSD")

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

convertToDfWithCols <- function(str,recSep,labelSep,valRm="",cols){
  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
  # vals <- sapply(vals, function(x)gsub(x,";",""))
  l <- as.list(vals)
  names(l) <- nms 
  df <- as.data.frame(l)
  cols.notfound <- cols[-which(cols %in% nms)]
  for ( i in cols.notfound){
    df[[i ]] <- NA
  }
  df
}


convertToList <- function(str,recSep,labelSep,valRm=""){
  records <- as.character(unlist(strsplit(str, recSep)))
  nms <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[1])))
  vals <- as.character(unlist(sapply(records, function(x)unlist(strsplit(x, labelSep))[2])))
  # vals <- sapply(vals, function(x)gsub(x,";",""))
  l <- as.list(vals)
  names(l) <- nms 
  l}

shortenIdVec <- function(gene_id, sep = "\\."){
  as.character(sapply(gene_id, function(x)unlist(strsplit(x,split=sep))[1]))
}



saveTranFibro <- function(fibroFile = fibroFile ){
  expr.df <-  getGeneExprForColMatch(query="Cells - Transformed fibroblasts",colName="SMTSD")
  exportAsTable(file=fibroFile,df=expr.df)
}



### DONT USE, really slow on 160k -> use cat/sed/awk scheme instead...
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
mergeENCODEexpr <- function (cyt.file=IMR.cyt.pap.tab,nuc.file=IMR.nuc.pap.tab) {
  cyt.df <- readInENCODEGtfToTab(file=cyt.file)
  nuc.df <- readInENCODEGtfToTab(file=nuc.file)
  cyt.cols <- colnames(cyt.df)
  nuc.cols <- colnames(nuc.df)
  colnames(cyt.df)[7:10] <- paste0(colnames(cyt.df)[7:10],".cyt")
  colnames(nuc.df)[7:10] <- paste0(colnamess(nuc.df)[7:10],".nuc")
  merge(cyt.df[6:10],nuc.df[5:10],by="transcript_id",all=TRUE)
}

loadInCytNucdata <- function(cyt.file=IMR.cyt.pap.tab,nuc.file=IMR.nuc.pap.tab,outfile=IMR.comb.expr ){
  imr.df <- mergeENCODEexpr(cyt.file, nuc.file)
  imr.idr.df <- imr.df[which(imr.df$IDR.nuc < 0.1 & imr.df$IDR.cyt < 0.1 &
                               !is.na(imr.df$IDR.cyt) & !is.na(imr.df$IDR.cyt) &
                               imr.df$IDR.nuc != 0 & imr.df$IDR.cyt != 0),]
  
  imr.idr.df$COMB.sum <- imr.idr.df$COMB.nuc + imr.idr.df$COMB.cyt
  getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"transcript_id"]),function(x)x[1]), ~ .id )
  lnc.expr.maxTransGeneCol.d <- suppressWarnings(ldply(split(imr.idr.df,imr.idr.df$gene_id),getTranscriptForMaxCols))
  trans.vec <- as.character(levels(lnc.expr.maxTransGeneCol.d$COMB.sum))
  imr.max.trans <- imr.idr.df[which(imr.idr.df$transcript_id %in% trans.vec),]
  exportAsTable(df =imr.max.trans,file=outfile)
}

getCombinedCytNucData <- function(file=IMR.comb.expr){
  df <- readInTable(file)
  lnc.genes <- getLncGenesV12()
  pc.genes  <- getPcGenesV12()
  lncIMR.df <- imr.df[which(imr.df$gene_id %in% lnc.genes$gene_id),]
  lncIMR.df$ratio <- lncIMR.df$COMB.cyt/lncIMR.df$COMB.nuc 
}




fetchEncodeDccFilesText <- function(filesTxt="~/data/wgEncodeCshlLongRnaSeqFiles.txt",
                                    filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab"){
  if (!file.exists(filesTxt)){
    download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/files.txt",dest=fileTxt)
  }
  
  d <- readLines(filesTxt)
  d.split <- lapply(d, function(x)strsplit(x,split="\t"))
  d.names <- sapply(d.split, function(x)unlist(x)[1])
  d.info <- sapply(d.split, function(x)unlist(x)[2])
  d.list <- lapply(d.info,function(x)convertToList(x,"; ", "=",";"))
  colAnnot <-  unique(do.call("c",lapply(d.list, function(x)names(x))))
  df <- do.call(rbind,lapply(d.info,function(x)convertToDfWithCols(x,"; ", "=",";",cols=colAnnot)))
  df$filename <- d.names
  
  exportAsTable(df=df,file=filesTxtTab)
}


getFilesTxt <- function(filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab",rnaExtract="longPolyA",gencodeVersion){
  if(!file.exists(filesTxtTab)){
    fetchEncodeDccFilesText()
  }
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df$localfile <- sapply(df$file, function(x)gsub(x,pattern="\\.gz", replacement=""))
  df$tabfile <- sapply(df$localfile, gtfToTabFilename)
  df$rnaExtractShort <- getRnaExtractShort(df$rnaExtract)
  df$cellShort <- sapply(df$cell, function(x)gsub(x, pattern="[-_]", replacement=""))
  df$headerMsg <- apply(df, 1,function(x)paste(x[["cellShort"]], x[["rnaExtractShort"]], x[["localization"]],sep="."))
  df$gencode <- tolower(str_extract(df$view,"(V[0-9][0-9]*)|DeNovo"))
  df$mappingTarget <- as.character(str_match(tolower(df$view),"gene|transcript|exon|fastq"))
  
  if(!missing(gencodeVersion)){
    df.out <- df[which(!is.na(df$tableName) & df$type == "gtf" & df$rnaExtract == rnaExtract & df$gencode == gencodeVersion),]
  } else {
    df.out <- df[which(!is.na(df$tableName) & df$type == "gtf" & df$rnaExtract == rnaExtract  ),]
    
  }
  
  df.out
}


getRnaExtractShort <- function(rnaExtractFn){
  vec <- ifelse(rnaExtractFn == "longPolyA", "lpa", rnaExtractFn)
  ifelse(vec == "longNonPolyA", "lnpa",vec)
  
}


getNucCytoExpr <- function(outfile="~/data/nucCytoEncodeExprFiles.tab",rnaExtract="longPolyA",gencodeVersion){
  if(!(rnaExtract %in% c("longPolyA", "longNonPolyA"))){
    stop("rnaExtract must be either <longPolyA> or <longNonPolyA>")
  }
  
  if(!missing(gencodeVersion)){
    if(!(gencodeVersion %in% c("v7", "v10"))){
      stop("rnaExtract must be either <v7> or <v10>")
    }
    df.gtf <- getFilesTxt(rnaExtract=rnaExtract,gencodeVersion=gencodeVersion)
  } else {
    df.gtf <- getFilesTxt(rnaExtract=rnaExtract)
    
    
  }
  nuc.df <- unique(df.gtf[which(df.gtf$localization == "nucleus"),])
  cyt.df <- unique(df.gtf[which(df.gtf$localization == "cytosol"), ])
  cytNuc.celltypes <- unique(cyt.df$cell[cyt.df$cell %in% nuc.df$cell])
  cytnuc.df <- df.gtf[which(!is.na(df.gtf$tableName) & df.gtf$type == "gtf" & df.gtf$rnaExtract == rnaExtract & (df.gtf$cell %in% cytNuc.celltypes)),]
  cytnuc.df
}


downloadCytNuc <- function(df){
  if(missing(df)){
    cytnuc.df <- df
  } else {
    cytnuc.df <- getNucCytoExpr()
  }
  gatherGzipFilesFromWeb(datadir=getDataFile(""),url.base= cshl.rnaseq.dir, url.file=cytnuc.df$filename)
}


gatherGzipFilesFromWeb <- function(datadir = local.datadir,
                                   url.base = cshl.rnaseq.dir,
                                   url.files){
  if (missing(datadir)){
    datadir <- getDataFile("")
  }
  if (!file.exists(datadir)){
    dir.create(datadir,recursive=TRUE)
  }
  
  print(datadir)
  url.vec <- paste0(url.base, url.files)
  temp.vec <- sapply(seq_along(url.vec),function(x) tempfile())
  
  invisible(sapply(seq_along(url.vec),function(x)download.file(url=url.vec[x],destfile=temp.vec[x])))
  
  invisible(sapply(seq_along(url.vec), 
                   function(i)write(readLines(con=gzfile(temp.vec[i]),n=-1),
                                    file=paste(datadir,"/",sub(x=url.files[i],pattern=".gz",replacement=""),sep=""))))
  
  
  
  invisible(sapply(seq_along(url.vec), function(i)unlink(temp.vec[i])))
}

gtfToTabFilename <- function(filename, datadir = local.datadir){
  nogz <- gsub(x=filename, pattern="\\.gz$",replacement="")
  tabfile <- gsub(x=nogz, pattern="\\.gtf$",replacement="\\.tab")
  
}

parseENCODEGtf <- function(gtfFile,tabFile,headerMsg){
  if(!file.exists(gtfFile)){
    stop("gtfFile not found")
  }
  if(file.exists(tabFile)){
    unlink(tabFile) 
  }
  header <- c("COMB", "RPKM1", "RPKM2", "IDR")
  header.title <- paste0(header,".", headerMsg)
  header.title.out <- do.call(paste,as.list(c("id", header.title)))
  tmp.cat <- tempfile()
  tmp.tab <- tempfile()
  system(paste("echo '", header.title.out,"' > ", tmp.cat, sep=""))
  system(paste("sed 's/[;\"]//g' ",gtfFile," | awk -F' ' '{print $12,$6,$14,$16,$18}' | sort -r> ",tmp.tab,sep=""))
  
  #cmd <- paste("cat <(echo ",header.title.out,") <(sed 's/[;\"]//g' ",gtfFile," | awk -F' ' '{print $12,$6,$14,$16,$18}' | sort) > ",tabFile,sep="")
  print(cmd)
  cmd <- paste("cat ",tmp.cat,tmp.tab,">",tabFile,sep=" ")
  system(cmd )
  unlink(tmp.cat)
  unlink(tmp.tab)
  #system(paste("cat ", tmp.cat, tmp.tab, ""))
}


systemJoin <- function(file1, file2, outfile){
  if(!file.exists("/usr/bin/join")){
    stop("cannot find join command in /usr/bin/join")
  }
  cmd <- paste("/usr/bin/join", file1, file2, ">",outfile,sep=" ")
  print(cmd)
  system(cmd)
}

removeCellTypes <- function(df, cellTypes){
  df[which(!(df$cellShort %in% cellTypes)),]
}
removeFileName <- function(df, filenames=cshl.remove.rnaseq){
  df[which(!(df$filename %in% cshl.remove.rnaseq)),]
}


mergeExprFiles <- function(outfile,gencodeVersion){
  if(missing(df)){
    df <- rbind(getNucCytoExpr(rnaExtract="longPolyA",gencodeVersion=gencodeVersion),
                getNucCytoExpr(rnaExtract="longNonPolyA",gencodeVersion=gencodeVersion))
    df <- removeFileName(df, filenames=cshl.remove.rnaseq)
    downloadCytNuc(df)
  }
  if(missing(outfile)){
    outfile <- getDataFile("encodeRnaSeqCshllpa.space")
  }
  
  df$localfile <- sapply(df$localfile, getDataFile)
  df$tabfile <- sapply(df$tabfile, getDataFile)
  
  #make short tab files
  for(i in seq_along(df$localfile)){
    parseENCODEGtf(gtfFile = df$localfile[i], tabFile = df$tabfile[i],headerMsg=df$headerMsg[i])
  }
  #dd <- read.csv(file=tabFile, sep=" ", stringsAsFactors=FALSE)
  
  tmp <- tempfile()
  tmp1 <- tempfile()
  file.copy(from = df$tabfile[1], to = tmp)
  for( i in 2:length(df$tabfile)){
    systemJoin(file1 = tmp, file2 = df$tabfile[i], outfile = tmp1)
    file.copy(from = tmp1, to = tmp)
  }
  file.copy(from=tmp, to = outfile)
  
  unlink(tmp)
  unlink(tmp1)
  
}







