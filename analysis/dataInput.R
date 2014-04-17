homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }

source(getFullPath("analysis/gencodeInput.R"))
source(getFullPath("analysis/clusterExecute.R"))


local.datadir <<- "/home/wespisea/data/"

cshl.remove.rnaseq <<- c("wgEncodeCshlLongRnaSeqNhekCellPapTranscriptGencV7Rep5.gtf.gz","wgEncodeCshlLongRnaSeqNhekCellPamTranscriptGencV7Rep5.gtf.gz")

rnaexpr.comb.v7 <<-  getDataFile("cshl-rnaSeq-v7-allcombined.space")
rnaexpr.comb.v7.cast <<-  getDataFile("cshl-rnaSeq-v7-allcombined-cast.space")
rnaexpr.comb.v7.melt <<-  getDataFile("cshl-rnaSeq-v7-allcombined-melt.space")
rnaexpr.comb.v7.maxTrans <<-  getDataFile("cshl-rnaSeq-v7-allcombined-maxTrans.space")


rnaexpr.comb.v10 <<-  getDataFile("cshl-rnaSeq-v10-allcombined.space")
rnaexpr.comb.v10.cast <<-  getDataFile("cshl-rnaSeq-v10-allcombined-cast.space")
rnaexpr.comb.v10.melt <<-  getDataFile("cshl-rnaSeq-v10-allcombined-melt.space")
rnaexpr.comb.v10.maxTrans <<-  getDataFile("cshl-rnaSeq-v10-allcombined-maxTrans.space")


cshl.rnaseq.dir <<- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"

gtex.annot <<- getDataFile("GTExSampleAnnot.tab")
gtex.gene.expr <<- getDataFile("GTExGeneRPKM.gct")
gtex.gene.head <<- getDataFile("GTExGeneRPKM.gct")
k562.nuc.pap <<- getDataFile("wgEncodeCshlLongRnaSeqK562NucleusPapTranscriptGencV7.gtf")
k562.nuc.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqK562NucleusPapTranscriptGencV7.tab")
k562.cyt.pap <<- getDataFile("wgEncodeCshlLongRnaSeqK562CytosolPapTranscriptGencV7.gtf")
k562.cyt.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqK562CytosolPapTranscriptGencV7.tab")
IMR.cyt.pap   <<-     getDataFile("wgEncodeCshlLongRnaSeqImr90CytosolPapTranscriptGencV10.gtf")
IMR.cyt.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqImr90CytosolPapTranscriptGencV10.tab")
IMR.nuc.pap <<-     getDataFile("wgEncodeCshlLongRnaSeqImr90NucleusPapTranscriptGencV10.gtf")
IMR.nuc.pap.tab <<- getDataFile("wgEncodeCshlLongRnaSeqImr90NucleusPapTranscriptGencV10.tab")
IMR.comb.expr   <<- getDataFile("ImrCytNucTopCombRnaExprTransByGene.tab")

fibroFile <<- getDataFile("transformedFibroblastsGTEx.tab")


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



getFilesTxtCytNucPapPam <- function(filesTxtTab="~/data/wgEncodeCshlLongRnaSeqFiles.tab"){
  if(!file.exists(filesTxtTab)){
    fetchEncodeDccFilesText()
  }
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))

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
  cytnuc.df <- df
  if(missing(df)){
    cytnuc.df <- getNucCytoExpr()
  }
  gatherGzipFilesFromWeb(datadir=getDataFile(""),url.base= cshl.rnaseq.dir, url.files=cytnuc.df$filename)
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
  ### TODO 
  # pull out seq_along(url.vec), and only download files that don't match checksums of localfiles
  # also, write force override
  download.index.all <- seq_along(url.vec)
  target.files <- sapply(seq_along(url.vec),function(i)paste(datadir,"/",sub(x=url.files[i],pattern=".gz",replacement=""),sep=""))
  file.exists <- sapply(seq_along(url.vec),function(i)file.exists(target.files[i]))
  download.index <- download.index.all
  download.index <- seq_along(url.vec)[which(file.exists == FALSE)]              
                         
  temp.vec <- sapply(download.index,function(x) tempfile())
  
  invisible(sapply(download.index,function(x)download.file(url=url.vec[x],destfile=temp.vec[x])))
  
  invisible(sapply(download.index, 
                   function(i)write(readLines(con=gzfile(temp.vec[i]),n=-1),
                                    file=target.files[i])))
  
  
  
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
 # print(cmd)
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

    df <- rbind(getNucCytoExpr(rnaExtract="longPolyA",gencodeVersion=gencodeVersion),
                getNucCytoExpr(rnaExtract="longNonPolyA",gencodeVersion=gencodeVersion))
    
    
    
    
  #  df <- removeFileName(df, filenames=cshl.remove.rnaseq)
  #  df <- removeCellTypes(df, cellTypes = "H1hESC")
  #  df <- removeCellTypes(df, cellTypes = "NHEK")
    
    downloadCytNuc(df=df)
  
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
  file.copy(from = df$tabfile[1], to = tmp, overwrite=TRUE)
  df$skip =FALSE
  for( i in 2:length(df$tabfile)){
    if(!checkCols(df$tabfile[i])){
      df$skip[i] = TRUE
      next
    }
    systemJoin(file1 = tmp, file2 = df$tabfile[i], outfile = tmp1)
    file.copy(from = tmp1, to = tmp, overwrite=TRUE)
  }
  file.copy(from=tmp, to = outfile, overwrite=TRUE)
  
  unlink(tmp)
  unlink(tmp1)
  checkAllCombined(allfile=outfile,df)
}


rep1 <- "wgEncodeCshlLongRnaSeqNhekCellPapTranscriptGencV7.gtf "
rep2 <- "wgEncodeCshlLongRnaSeqNhekCellPapTranscriptGencV7Rep5.gtf"

naVecToVal <- function(vec,val){
  vec[is.na(vec)] <- val
  vec
}

checkCols <- function (allfile) {
  colpipe <- pipe(paste("cat",allfile,"| awk -F ' ' '{print NF}'"))
  colCounts <- as.numeric(readLines(colpipe))
  close(colpipe)
  if(!sum(colCounts == colCounts[1]) == length(colCounts))
    return(FALSE)
  return(TRUE)
}


checkAllCombined <- function(allfile = rnaexpr.comb.v7,df){
  pass = TRUE
  all.df <- read.csv(allfile,stringsAsFactors=FALSE,sep=" ", nrows=10000)
  checkCols(allfile)
  
  result <- list()
  for(line in seq_len(nrow(df))){
    if(!identical(df$skip[line], TRUE)){
    for(reading in c("RPKM1", "RPKM2", "COMB", "IDR")){
      testLine <- df[line,]
      comb.vec <- all.df[[paste(reading,testLine$headerMsg,sep=".")]]
      gtf.df <- read.csv(file=testLine$tabfile, sep = " ", stringsAsFactors=FALSE, nrows=10000)
      col <- paste(reading,testLine$headerMsg,sep=".")
      tabvalues <- naVecToVal(gtf.df[[col]],-1)
      combvalues <- naVecToVal(all.df[[col]],-1)
      result[[col]] <- (sum(tabvalues == combvalues ) == length(tabvalues))
      }
    }
  }
  if(any(as.logical(result) == FALSE)){
    return(FALSE)
  }else {
    return(TRUE)
  }
  
}

allExprCombList <- function(names.df){
   l <- list(comb = names.df[grep("COMB",names.df)],
       RPKM1 = names.df[grep("RPKM1",names.df)],
       RPKM2 = names.df[grep("RPKM2",names.df)],
       IDR = names.df[grep("IDR",names.df)],
       lpa = names.df[grep("lpa",names.df)],
       lnpa = names.df[grep("lnpa",names.df)],
       cell = names.df[grep("cell",names.df)],
       cytosol = names.df[grep("cytosol",names.df)],
       nucleus = names.df[grep("nucleus",names.df)])
   
   celltypes <- unique(as.character(sapply(names.df[-1], function(x)unlist(strsplit(x,"\\.")[[1]])[2])))
   cell.list <- lapply(celltypes,function(x)names.df[grep(x,names.df)])
   names(cell.list) <- celltypes
   out.list <- c(l,cell.list)
   out.list
}

combLabelVecPart <- function(labelVec, comp){
  component <- c("measure", "cell", "pulldown", "localization")
  index <- which(component %in% comp)
  if(identical(index,integer(0))){
    stop("comp not in c(measure,cell,pulldown,localization)")
  }
  as.character(sapply(labelVec, function(label)unlist(strsplit(label,"\\.")[[1]])[index]))
}

castOnMeasure <- function(df = df.genes.melt){
  #df.skinny <- df[c("id","gene_id","value","measure")]
  d <- do.call(cbind, split(df[c("value","id")],df$measure))
  c.1 <- sum(d$COMB.id == d$RPKM1.id)/(dim(d)[1])
  c.2 <-  sum(d$RPKM2.id == d$RPKM1.id)/(dim(d)[1])
  c.3 <-  sum(d$RPKM2.id == d$IDR.id)/(dim(d)[1])
  if(identical(c.1,1) && identical(c.2,1) && identical(c.3,1)){
    d.out <- d[grep("value", names(d))] 
    colnames(d.out) <- gsub(colnames(d.out), pattern=".value",replacement="")
    d1 <- cbind(df[which(df$measure == (df$measure)[1]),], d.out)
    d1$value <- NULL
    d1$measure <- NULL
    return(d1)
  }
  
}


processExprComb <- function(gencodeVersion="v7"){
  file.expr <- get(paste("rnaexpr.comb", gencodeVersion,sep="."),envir=parent.frame())
  if (!file.exists(file.expr)){
    stop("cannot find expr file...aborating")
  }
  df <- read.csv(file=file.expr, sep = " ", stringsAsFactors=FALSE)
  col.list <- allExprCombList(names(df))
  pc.genes <- getGencodeAnnot(biotype="pc", gencodeVersion)
  lnc.genes <- getGencodeAnnot(biotype="lnc", gencodeVersion)
  df.pcLnc <- df[which(df$id %in% c(pc.genes$transcript_id, lnc.genes$transcript_id)),]
  df.genes <- merge(df.pcLnc, rbind(pc.genes[c("transcript_id", "gene_id")],lnc.genes[c("transcript_id", "gene_id")]),by.x="id",by.y="transcript_id",all.x=TRUE)
 
#  d1 <- df.genes
#  df.genes <- df.genes[1:1000,]
  #df.rpkmsum <- df.genes[c("id","gene_id",rpkmsum)]
  df.genes.melt <- melt(df.genes, id.vars=c("id", "gene_id"))
  df.genes.melt$var = levels(df.genes.melt$variable)[df.genes.melt$variable]
  df.genes.melt$measure <- combLabelVecPart(df.genes.melt$var,"measure")
  df.genes.melt$cell <- combLabelVecPart(df.genes.melt$var,"cell")
  df.genes.melt$pulldown <- combLabelVecPart(df.genes.melt$var,"pulldown")
  df.genes.melt$localization <- combLabelVecPart(df.genes.melt$var,"localization")
  
  df.genes.melt$variable <- NULL
  df.genes.melt$var <- NULL
  file.melt <- get(paste("rnaexpr.comb", gencodeVersion,"melt",sep="."),envir=parent.frame())
  exportAsTable(df.genes.melt,file.cast)
  
  
  
  df.genes.cast <- castOnMeasure(df.genes.melt)
  #df.genes.cast <- dcast(df.genes.melt, ... ~ measure)
  df.genes.cast$RPKMsum <- df.genes.cast$RPKM1 +  df.genes.cast$RPKM2
  file.cast <- get(paste("rnaexpr.comb", gencodeVersion,"cast",sep="."),envir=parent.frame())
  exportAsTable(df.genes.cast,file.cast)
  
  df.genes.cast$dd <- factor(paste0(df.genes.cast$gene_id, df.genes.cast$cell,df.genes.cast$pulldown,df.genes.cast$localization))
  dgc.table <- tbl_df(df.genes.cast)
  dgc.table.dd <- group_by(dgc.table, dd)
  df.genes.trans <- as.data.frame(filter(dgc.table.dd,RPKMsum == max(RPKMsum) ))
  
  #df.genes.trans <- ddply(df.genes.cast, .(dd),
   #                       function(x)x[which(x$RPKMsum == max(x$RPKMsum))[1],])
  df.genes.trans$dd <- NULL
  file.maxTrans <- get(paste("rnaexpr.comb", gencodeVersion,"maxTrans",sep="."),envir=parent.frame())
  exportAsTable(file.maxTrans,file.maxTrans)
  
  return(1)
  
  df.genes.trans$val <- paste0(df.genes.trans$cell,".",df.genes.trans$pulldown,".",df.genes.trans$localization)
  df.genes.trans$id <- NULL
  df.genes.trans$cell <- NULL
  df.genes.trans$pulldown <- NULL
  df.genes.trans$localization <- NULL
  df.genes.trans.melt <- melt(df.genes.trans, id.vars=c("val", "gene_id"))
  df.genes.trans.melt$var <- levels(df.genes.trans.melt$variable)[df.genes.trans.melt$variable]
  df.genes.trans.melt$variable <- NULL
  df.genes.trans.melt$colname <- paste0(df.genes.trans.melt$var,".",df.genes.trans.melt$val)
  df.genes.trans.melt$var <- NULL
  df.genes.trans.melt$val <- NULL
  df.genes.maxTrans <- dcast(df.genes.trans.melt, ... ~ colname)
  
  #getTranscriptForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)t1[which(x == max(x)),"id"]),function(x)x[1]), ~ .id )
  #getValForMaxCols <- function(t1)cast(ldply(apply(t1,2,function(x)x[which(x == max(x))]),function(x)x[1]), ~ .id )
  #lnc.expr.maxTransGeneCol.d <- suppressWarnings(ldply(split(imr.idr.df,imr.idr.df$gene_id),getValForMaxCols))
}


createRnaSeqMapFile <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  all.equal(df.comb$read2.localization,df.comb$read1.localization)
  all.equal(df.comb$read2.rnaExtract,df.comb$read1.rnaExtract)
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.sam")
  paste0("STAR  --readFilesCommand zcat --readFilesIn ", rnaseqdir,df.comb$read1.filename, " ",rnaseqdir,df.comb$read2.filename,
         " > ",rnaseqdir,df.comb$output,sep=" ")
}


createRnaSeqSegemehl <- function(){
  # ./segemehl.x -i chr1_2_3.idx -d chr1.fa chr2.fa chr3.fa -q myreads.fa --threads 8 > mymap.sam
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  all.equal(df.comb$read2.localization,df.comb$read1.localization)
  all.equal(df.comb$read2.rnaExtract,df.comb$read1.rnaExtract)
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  seg <-  "~/bin/segemehl_0_1_7/segemehl/segemehl.x "
  seg.idx <- "/project/umw_zhiping_weng/wespisea/rna-seq/GRCh37.p13.genome.idx"
  seg.fa <- "/project/umw_zhiping_weng/wespisea/rna-seq/GRCh37.p13.genome.fa"
  
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="seg.sam")
  df.comb$bam <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="seg.bam")
  
  cat(paste0(seg,"-t 16 -i ",seg.idx, " -d " ,seg.fa ," -q ", rnaseqdir,df.comb$read1.filename, " -p ",rnaseqdir,df.comb$read2.filename,
         " > ",rnaseqdir,df.comb$output,"\n",sep=" "),file="~/sandbox/segRun.sh")
  
  #paste0("samtools view -bS ",rnaseqdir,df.comb$output, " > ",rnaseqdir,df.comb$bam), file="~/sandbox/segRunBam.sh")

   write(generateSamToolsCoversionStar(df.comb$bare), file="~/sandbox/segRunBam.sh")
}

createRnaSeqSegemehlRest <- function(){
  # ./segemehl.x -i chr1_2_3.idx -d chr1.fa chr2.fa chr3.fa -q myreads.fa --threads 8 > mymap.sam
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell != "K562" & cell != "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  all.equal(df.comb$read2.localization,df.comb$read1.localization)
  all.equal(df.comb$read2.rnaExtract,df.comb$read1.rnaExtract)
   
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="seg.sam")
  df.comb$bam <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="seg.bam")
  
  cat(getSegCmdMap(df.comb$read1.filename,df.comb$read2.filename,df.comb$output),file="~/sandbox/segRunRest.sh")
  
  #paste0("samtools view -bS ",rnaseqdir,df.comb$output, " > ",rnaseqdir,df.comb$bam), file="~/sandbox/segRunBam.sh")
  write(generatePicardCoversion(df.comb$bare), file="~/sandbox/segRunRest.sh")
  generatePicardCoversion(df.comb$bare)
}




getSegCmdMap <- function(rd1,rd2,output){
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  seg <-  "~/bin/segemehl_0_1_7/segemehl/segemehl.x "
  seg.idx <- "/project/umw_zhiping_weng/wespisea/rna-seq/GRCh37.p13.genome.idx"
  seg.fa <- "/project/umw_zhiping_weng/wespisea/rna-seq/GRCh37.p13.genome.fa"
  paste0(seg,"-t 16 -i ",seg.idx, " -d " ,seg.fa ," -q ", rd1, " -p ",rd2,
         " > ",output,sep=" ")
  
}




gbToBytes <- function(gb){
  gb * (1024)^3
}
generateSamToolsCoversion <- function(filebase,dir=rnaseqdir,tag="seg"){
  size <- gbToBytes(40)
  cmd <- "~/bin/samtools view -bS test.seg.sam > test.seg.bam;;~/bin/samtools sort -m 42949672960 test.seg.bam test.seg.sort;;~/bin/samtools index test.seg.sort.bam test.seg.bai;;~/bin/samtools idxstats test.seg.bai > test.seg.bai.stats"
  as.character(unlist(sapply(filebase, function(filename)gsub(x=gsub(x=cmd,pattern="test", replacement=file.path(dir,filename)),pattern="seg", replacement=tag))))
}

#.star.samAligned.out.sam

generateSamToolsCoversionStar <- function(filebase,dir=rnaseqdir,tag="seg"){
  size <- gbToBytes(40)
  cmd <- "~/bin/samtools view -bS test.star.samAligned.out.sam > test.seg.bam;;~/bin/samtools sort -m 42949672960 test.seg.bam test.seg.sort;;~/bin/samtools index test.seg.sort.bam test.seg.bai;;~/bin/samtools idxstats test.seg.bai > test.seg.bai.stats"
  as.character(unlist(sapply(filebase, function(filename)gsub(x=gsub(x=cmd,pattern="test", replacement=file.path(dir,filename)),pattern="seg", replacement=tag))))
}

generatePicardCoversionStar <- function(filebase,dir=rnaseqdir,tag="seg"){
  size <- gbToBytes(40)
  cmd <- "java -jar -Xmx30g /share/pkg/picard/1.96/SamFormatConverter.jar INPUT=test.star.samAligned.out.sam OUTPUT=test.seg.bam;;java -jar -Xmx30g  /share/pkg/picard/1.96/SortSam.jar SORT_ORDER=coordinate INPUT=test.seg.bam OUTPUT=test.seg.sort.bam;;java -jar -Xmx30g  /share/pkg/picard/1.96/BuildBamIndex.jar INPUT=test.seg.sort.bam OUTPUT=test.seg.bai;;java -jar -Xmx30g  /share/pkg/picard/1.96/BamIndexStats.jar INPUT=test.seg.bai OUTPUT=test.seg.bai.stats"
  as.character(unlist(sapply(filebase, function(filename)gsub(x=gsub(x=cmd,pattern="test", replacement=file.path(dir,filename)),pattern="seg", replacement=tag))))
}
generatePicardCoversion <- function(filebase,dir=rnaseqdir,tag="seg"){
  size <- gbToBytes(40)
  cmd <- "java -jar -Xmx30g /share/pkg/picard/1.96/SamFormatConverter.jar INPUT=test.seg.sam OUTPUT=test.seg.bam;;java -jar -Xmx30g  /share/pkg/picard/1.96/SortSam.jar SORT_ORDER=coordinate INPUT=test.seg.bam OUTPUT=test.seg.sort.bam;;java -jar -Xmx30g  /share/pkg/picard/1.96/BuildBamIndex.jar INPUT=test.seg.sort.bam OUTPUT=test.seg.bai;;java -jar -Xmx30g  /share/pkg/picard/1.96/BamIndexStats.jar INPUT=test.seg.bai OUTPUT=test.seg.bai.stats"
  as.character(unlist(sapply(filebase, function(filename)gsub(x=gsub(x=cmd,pattern="test", replacement=file.path(dir,filename)),pattern="seg", replacement=tag))))
}



getRead12ForK562andGM <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.sam")
  df.comb$bam <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.bam")
  
  genStarAlignCmd(df.comb$read1.filename, df.comb$read2.filename,df.comb$output)
}

getRead12ForK562andGM.spikeIn <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.sam")
  df.comb$bam <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.bam")
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  o <- paste0(genStarAlignCmdSpikeIN(df.comb$read1.filename, df.comb$read2.filename,df.comb$output),
  "; ", generateSamToolsCoversion(df.comb$bare,tag="star"))
  write(o, file="~/sandbox/starSpikeBam.sh")
  
  o.bam <-  generatePicardCoversionStar(df.comb$bare,tag="star")
  write(o.bam, file="~/sandbox/starSpikeBam.sh")
  
}


generateRestOfStarSpikeIN <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell != "K562" | cell != "GM12878"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$output <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.sam")
  df.comb$bam <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern="fastq.gz",replacement="star.bam")
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  o <- paste0(genStarAlignCmdSpikeIN(df.comb$read1.filename, df.comb$read2.filename,df.comb$output),
              "; ", generatePicardCoversionStar(df.comb$bare,tag="star"))
  write(o, file="~/sandbox/starSpikeBamRest.sh")
  
}


starGenerateGenomeCommand <- function(){
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDir/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  #STAR --runMode genomeGenerate --genomeDir genomepath --genomeFastaFiles  genomepath/genome.fa  
  # --sjdbGTFfile genomepath/genes.gtf --sjdbOverhang 100 --runThreadN 8
  paste0("STAR --runMode genomeGenerate --genomeDir ", starGenomeDir, " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf, " --sjdbOverhand 75 --runThreadN 16") 
}
starGenerateGenomeCommand.spike <- function(){
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDirSpikeIn/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.spikeIn.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  #STAR --runMode genomeGenerate --genomeDir genomepath --genomeFastaFiles  genomepath/genome.fa  
  # --sjdbGTFfile genomepath/genes.gtf --sjdbOverhang 100 --runThreadN 8
  star.cmd <- paste0("STAR --runMode genomeGenerate --genomeLoad LoadAndRemove --genomeDir ", starGenomeDir, " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf, " --sjdbOverhang 75 --runThreadN 16") 
  cat("perl ~/bin/runJob.pl -c 20 -m 40960 -W 600 -Q short -i \"",star.cmd,"\"",sep="" )

}

# 
genStarAlignCmdSpikeIN <- function(rd1,rd2,outfile){
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDirSpikeIn/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.spikeIn.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  
  paste0("STAR --runMode alignReads --genomeLoad LoadAndRemove --readFilesCommand zcat --runThreadN 16", 
         " --readFilesIn ", paste0(rnaseqdir, rd1), " ", paste0(rnaseqdir, rd2), 
         " --genomeDir ", starGenomeDir,
         " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf,
         " --outFileNamePrefix ", paste0(rnaseqdir, outfile)) 
}


genStarAlignCmdSpikeIN.paramFile <- function(rd1,rd2,outfile){
  paramFile <- "/home/aw30w/log/params/parametersENCODElong_AWmod.txt"
  paste0("STAR --runMode alignReads ", 
         " --readFilesIn ", rd1, " ", rd2, 
         " --outFileNamePrefix ", outfile,
         " --parametersFiles ", paramFile)
}



genStarAlignCmd <- function(rd1,rd2,outfile){
  rnaseqdir <- "/project/umw_zhiping_weng/wespisea/rna-seq/"
  starGenomeDir <- paste(rnaseqdir,"starGenomeDir/",sep="")
  genomeFasta <- paste(rnaseqdir,"GRCh37.p13.genome.fa",sep="")
  annotationGtf <- paste(rnaseqdir,"gencode.v19.annotation.gtf",sep="")
  
  paste0("STAR --runMode alignReads --genomeLoad LoadAndRemove --runThreadN 16", 
         " --readFilesIn ", paste0(rnaseqdir, rd1), " ", paste0(rnaseqdir, rd2), 
         " --genomeDir ", starGenomeDir,
         " --genomeFastaFiles ", genomeFasta,
         " --sjdbGTFfile ", annotationGtf,
         " --outFileNamePrefix ", paste0(rnaseqdir, outfile)) 
}



downloadBamBaiFile <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.cytNuc.restFastq <- subset(df, type == "fastq" & (localization == "nucleus" | localization == "cytosol") & cell != "K562" & cell != "GM12878")
  df.cytNuc.bam <- subset(df, type == "bam" & (localization == "nucleus" | localization == "cytosol") & cell != "K562" & cell != "GM12878")
  filenames <- c(df.cytNuc.restFastq$filename, df.cytNuc.bam$filename)
  remote.site <- cshl.rnaseq.dir
  rna.remote <- file.path("/project/umw_zhiping_weng/wespisea/","rna-seq/")
  web.file <- paste0(remote.site,filenames)
  #sapply(web.file, function(w)url.exists(w))
  hpc.download(url.vec=paste0(remote.site,filenames), target.vec=paste0(rna.remote,filenames))
  hpc.download.seq.cont(url.vec=paste0(remote.site,filenames), target.vec=paste0(rna.remote,filenames))
}


getSpikeIns <- function(){
  spikeIn <- read.fasta(file=url("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/supplemental/wgEncodeCshlLongSpikeins.fasta"))
  lenList <- sapply(spikeIn, function(x)length(x))
  nms <- names(lenList)
  val <- as.numeric(unlist(lenList))
  valp1 <- val + 1
  nms.trans <- paste0(nms,"_trans")
  spikeInbed12 <- paste0(paste0(nms,"\t",1,"\t",valp1,"\t",nms.trans,"\t",0,"\t","+","\t",1,"\t",valp1,"\t",0,"\t",1,"\t",val,",\t","0,"),"\n")
  cat(spikeInbed12, file="~/sandbox/spikeIn14.bed12")
}

getSpikeInsGTF <- function(){
  spikeIn <- read.fasta(file=url("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/supplemental/wgEncodeCshlLongSpikeins.fasta"))
  lenList <- sapply(spikeIn, function(x)length(x))
  nms <- names(lenList)
  val <- as.numeric(unlist(lenList))
  valp1 <- val + 1
  nms.trans <- paste0(nms,"_trans")
  nms.gene  <- paste0(nms,"_gene")
  nms.exon  <- paste0(nms,"_exon")
  
  spikeInGTF.gene <- paste0(nms,"\t","NISTSpikeIn","\t","gene","\t",1,"\t",valp1,"\t",".","\t","+","\t",".","\t",
                              "gene_id \"",nms.gene,"\";"," transcript_id \"",nms.trans,"\";")
  spikeInGTF.trans <- paste0(nms,"\t","NISTSpikeIn","\t","transcript","\t",1,"\t",valp1,"\t",".","\t","+","\t",".","\t",
                            "gene_id \"",nms.gene,"\";"," transcript_id \"",nms.trans,"\";")
  spikeInGTF.exon <- paste0(nms,"\t","NISTSpikeIn","\t","exon","\t",1,"\t",valp1,"\t",".","\t","+","\t",".","\t",
                             "gene_id \"", nms.gene, "\";"," transcript_id \"",nms.trans,"\"; ", " exon_id \"",nms.exon,"\";")
  gtfLines <- paste0(spikeInGTF.gene,"\n",spikeInGTF.trans,"\n",spikeInGTF.exon,"\n")
  
  
  
  cat(gtfLines, file="~/sandbox/spikeIn14.gtf")
  # perl -pi -e 's/^ //g'  ~/sandbox/spikeIn14.gtf 
  # scp ~/sandbox/spikeIn14.gtf aw30w@ghpcc06.umassrc.org:/project/umw_zhiping_weng/wespisea/rna-seq/
  # cat /project/umw_zhiping_weng/wespisea/rna-seq/gencode.v19.annotation.gtf /project/umw_zhiping_weng/wespisea/rna-seq/spikeIn14.gtf > /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.NIST14SpikeIn.gtf

  
  
  
}



getBedopsIntersect <- function(infile,tag){
  
if(tag == "star"){

  # samtools view -bS test.star.samAligned.out.sam > test.star.bam;;bedtools bamtobed -i test.star.bam > test.star.bed
  genbed <- "samtools view -bS test.star.samAligned.out.sam > test.star.bam;;bedtools bamtobed -i test.star.bam > test.star.bed"
  sortbed <- "/home/aw30w/bin/sort-parallel.sh test.star.bed test.star_sort.bed"
  exon<-"/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exon.bed > test.exon.found"
  intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intron.bed > test.intron.found"
  gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
  cmd <- paste0(c(genbed,sortbed,exon,intron,gene),collapse=";;")
} else {
  genbed <-"/home/aw30w/bin/bedops-2.4.1/bin/sam2bed < test.seg.sam > test.seg.bed"
  sortbed <- "/home/aw30w/bin/sort-parallel.sh test.seg.bed > test.seg_sort.bed"
  exon <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exon.bed > test.exon.found"
  intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intron.bed > test.intron.found"
  gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
  cmd <- paste0(c(genbed,sortbed,exon,intron,gene),collapse=";;")
  
}
as.character(unlist(sapply(infile, function(xx)gsub(x=cmd,pattern="test",replacement=xx))))

}

getBedopsIntersectStartWithBed <- function(infile,tag){
  
  if(tag == "star"){
    
    # samtools view -bS test.star.samAligned.out.sam > test.star.bam;;bedtools bamtobed -i test.star.bam > test.star.bed
    sortbed <- "/home/aw30w/bin/sort-parallel.sh test.star.bed  test.star_sort.bed"
    exon<-"/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exon.bed > test.exon.found"
    intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intron.bed > test.intron.found"
    gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
    cmd <- paste0(c(sortbed,exon,intron,gene),collapse=";;")
  } else {
    sortbed <- "/home/aw30w/bin/sort-parallel.sh test.seg.bed test.seg_sort.bed"
    exon <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exon.bed > test.exon.found"
    intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intron.bed > test.intron.found"
    gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
    cmd <- paste0(c(sortbed,exon,intron,gene),collapse=";;")
    
  }
  as.character(unlist(sapply(infile, function(xx)gsub(x=cmd,pattern="test",replacement=xx))))
  
}

getBedopsIntersect <- function(infile,tag){
  
  if(tag == "star"){
    
    # samtools view -bS test.star.samAligned.out.sam > test.star.bam;;bedtools bamtobed -i test.star.bam > test.star.bed
  #  cmd <- "/home/aw30w/bin/sort-parallel.sh test.star.bed  test.star_sort.bed"
    exon <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exonMerge.bed > test.exon.found"
    intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intronExDiff.bed > test.intron.found"
    gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.star_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
    cmd <- paste0(c(exon,intron,gene),collapse=";;")
    
    
    
  } else {
    #sortbed <- "/home/aw30w/bin/sort-parallel.sh test.seg.bed test.seg_sort.bed"
    # ../gencodeV19.exonMerge.bed > ../gencodeV19.intronExDiff.bed
    
    exon <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.exonMerge.bed > test.exon.found"
    intron <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.intronExDiff.bed > test.intron.found"
    gene <- "/home/aw30w/bin/bedops-2.4.1/bin/bedops --element-of -50%  test.seg_sort.bed /project/umw_zhiping_weng/wespisea/rna-seq/gencodeV19.bed > test.gene.found"
    cmd <- paste0(c(exon,intron,gene),collapse=";;")
    
  }
  as.character(unlist(sapply(infile, function(xx)gsub(x=cmd,pattern="test",replacement=xx))))
  
}

generateStarBedops<- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  o <- paste0(getBedopsIntersect(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),tag="star"))
  write(o, file="~/sandbox/starToBed.sh")
  # cat ~/bin/starToBed.sh | xargs -I{}  perl ~/bin/runJob.pl -c 2 -m 20480 -W 600 -Q short -t "star2Bed" -i "{}"
}


sortStarBed <- function(base){
  cmd1 <- "/home/aw30w/bin/bedops-2.4.1/bin/sort-bed --max-mem 40G --tmpdir /project/umw_zhiping_weng/wespisea/tmp test.star.bed > test.star_sort.bed"
  #cmd2 <- "/home/aw30w/bin/bedops-2.4.1/bin/sort-bed --max-mem 20G --tmpdir /project/umw_zhiping_weng/wespisea/tmp test.star.bed test.star_sorted.bed"
  #cmd3 <- "/home/aw30w/bin/bedops-2.4.1/bin/sort-bed --max-mem 20G --tmpdir /project/umw_zhiping_weng/wespisea/tmp test.star.bed test.star_sorted.bed"
  as.character(unlist(sapply(base, function(xx)gsub(x=cmd1,pattern="test",replacement=xx))))
}

generateStarBedops<- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol") & (cell == "K562" | cell == "GM12878"))
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  
  o <- genStarAlignCmdSpikeIN.paramFile(rd1=file.path(rnaseqdir,df.comb$read1.filename),
                                        rd2=file.path(rnaseqdir,df.comb$read2.filename),
                                        outfile=file.path(rnaseqdir,"starSpikeIn",df.comb$bare))
  write(o, file="~/sandbox/runStar.sh")
  scpFile(file.local="~/sandbox/runStar.sh", dir.remote="~/bin/")
  # cat ~/bin/runStar.sh | xargs -I{}  perl ~/bin/runJob.pl -c 16 -m 3072 -W 600 -Q short -t "runStar" -i "{}"
  df.comb$starAln <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star.samAligned.out.sam"))
  
  
  cmd1 <- "samtools view -bS test.star.samAligned.out.sam -o test.star.bam;;samtools sort -m 171798691840 test.star.bam test.star_sort"
  o1 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  write(o1, file="~/sandbox/starBamSamtools.sh")
  scpFile(file.local="~/sandbox/starBamSamtools.sh", dir.remote="~/bin/")
  df.comb$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam"))
# cat ~/bin/starBamSamtools.sh | perl -e 'while(<>){@c = split(/ +/,$_);$file = $c[scalar(@c)-1]; chomp($file); print $_ if !-e $file.".bam" }' | xargs -I{}  perl ~/bin/runJob.pl -c 2 -m 183840 -W 600 -Q short -t "samtools" -i "{}"  
  #cmd2 <- "java -jar -Xmx70g /share/pkg/picard/1.96/SamFormatConverter.jar INPUT=test.star.samAligned.out.sam OUTPUT=test.starPic.bam;;java -jar -Xmx30g  /share/pkg/picard/1.96/SortSam.jar SORT_ORDER=coordinate INPUT=test.starPic.bam OUTPUT=test.starPic.bam"
  #o2 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd2,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  #write(o2, file="~/sandbox/starBamPicard")
  #scpFile(file.local="~/sandbox/starBamPicard", dir.remote="~/bin/")
  #df.comb$picSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".starPic_bam"))
  
  # cat ~/bin/starBamPicard | xargs -I{}  perl ~/bin/runJob.pl -c 2 -m 40960 -W 600 -Q short -t "picard" -i "{}"
  
   cmd3 <- "samtools index test.star_sort.bam"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  write(o3, file="~/sandbox/stIdxStar")
  scpFile(file.local="~/sandbox/stIdxStar", dir.remote="~/bin/")
  df.comb$bamBai <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam.bai"))
  # cat ~/bin/stIdxStar | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 4096 -W 600 -Q short -t bamIndexTest -i "{}"
  
  
  
  cmd4 <- "samtools view -h test.star_sort.bam -o test.star_sort.sam"
  o4 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd4,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  write(o4, file="~/sandbox/sortBam2Sam")
  scpFile(file.local="~/sandbox/sortBam2Sam", dir.remote="~/bin/")
  df.comb$samSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.sam"))
  # catMissingLast.sh ~/bin/sortBam2Sam | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 8192 -W 600 -Q short -t sortBam2sam -i "{}"
  
  
  
  cmd5.1 <- "bedtools multicov -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.bed > COUNT_FILE.multicov.allTrans.uniqReads"
  cmd5.2 <- "bedtools multicov -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.singleTrans.bed > COUNT_FILE.multicov.singleTrans.uniqReads"
  cmd5.3 <- "bedtools multicov -D -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.bed > COUNT_FILE.multicov.allTrans.dupReads"
  cmd5.4 <- "bedtools multicov -D -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.singleTrans.bed > COUNT_FILE.multicov.singleTrans.dupReads"
  cmd5.5 <- "bedtools multicov  -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.grace.bed > COUNT_FILE.multicov.grace"

  o5 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste(cmd5.1,cmd5.2,cmd5.3,cmd5.4,cmd5.5,sep="\n"),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)),
                                                                        pattern="COUNT_FILE",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/multicovCounts/",filename)))))
  write(o5, file="~/sandbox/multicovCounts")
  scpFile(file.local="~/sandbox/multicovCounts", dir.remote="~/bin/")
  df.comb$multicovBare <- file.path(rnaseqdir,"starSpikeIn/multicovCounts/",df.comb$bare)
  # catMissingLast.sh ~/bin/multicovCounts | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 8192 -W 600 -Q short -t multicov -i "{}"
  # cat ~/bin/multicovCounts | grep "multicov.grace"| xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 8192 -W 600 -Q short -t multicov -i "{}"
  
  df.comb$multicovGrace <- paste0(rnaseqdir,"/starSpikeIn/multicovCounts/",df.comb$bare,".multicov.grace")
  df.tiny <- df.comb[c("read1.cell","read1.localization","read1.rnaExtract","read1.replicate","multicovGrace")]
  colnames(df.tiny) <- c("cell","localization","rnaExtract","replicate","multicovGrace")
  grace.df <- merge(df.tiny,expand.grid(region = c("cds", "noncoding", "utr3", "utr5","intron"), multicovGrace=df.comb$multicovGrace),by="multicovGrace")
  grace.df$cmdOut <- paste0(grace.df$multicovGrace,".",grace.df$region,".summary")
  grace.df$cmd <- paste0("cat ",grace.df$multicovGrace," | grep ",grace.df$region," | awk '{NR += \\$7;} END {print NR}'")
  
  counts <- c()
  for(i in seq_along(grace.df$cmd)){
    counts[i] <- as.numeric(hpc.system(grace.df$cmd[i]))
  }
  grace.df$counts <- counts
  exportAsTable(df=grace.df,file=getFullPath("data/multicovCountsGraceFeatures.tab"))
  
  
  
  cmd6.1 <- paste("cat /project/umw_zhiping_weng/wespisea/flux-capicitor/v19annotParamsNIST14 >> test.flux.params && echo STDOUT_FILE test.flux.out >> test.flux.params &&",
                  "echo STDERR_FILE test.flux.err >> test.flux.params && echo STATS_FILE test.flux.stats >> test.flux.params && echo COVERAGE_FILE test.flux.coverage ",
                  ">> test.flux.params")
  cmd6.2 <- paste("/home/aw30w/bin/flux-capacitor-1.6.1/bin/flux-capacitor --threads 24 -p test.flux.params",
                "-i  xxTESTINPUTxx.star_sort.bam -o test.flux.output")
  
  o6 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste0(cmd6.1,";;",cmd6.2),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14",filename)),
                                                                        pattern="xxTESTINPUTxx",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/",filename)  ))))
  write(o6, file="~/sandbox/fluxCap")
  scpFile(file.local="~/sandbox/fluxCap", dir.remote="~/bin/")
  # perl /home/aw30w/bin/runTask.pl -f ~/bin/fluxCap -c 24 -m 1024 -W 600 -Q short -t flux
  
  
  
  
  # /home/aw30w/bin/sjcount-master/sjcount -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file]
  cmd7 <- paste0("/home/aw30w/bin/sjcount-master/sjcount  -quiet  -bam  xxTESTINPUTxx.star_sort.bam -ssj test.ssj -ssc test.scc -log test.log" )
  o7 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste0(cmd7),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn/ssjCount",filename)),
                                                                        pattern="xxTESTINPUTxx",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/",filename)  ))))
  write(o7, file="~/sandbox/ssjcountRun")
  scpFile(file.local="~/sandbox/ssjcountRun", dir.remote="~/bin/")
  # perl /home/aw30w/bin/runTask.pl -f ~/bin/ssjcountRun -c 2 -m 8192 -W 600 -Q short -t ssj
  
  
  
  
  cmd5.1 <- "python reads-in-features.py  --sam=test.star_sort.sam --gff=/project/umw_zhiping_weng/wespisea/gtf/gencode.v19.annotation.pc.gtf  --label=mRNA --stranded=TRUE --outprefix=COUNT_FILE.mRNA."
  cmd5.2 <- "python reads-in-features.py  --sam=test.star_sort.sam --gff=/project/umw_zhiping_weng/wespisea/gtf/gencode.v19.long_noncoding_RNAs.gtf  --label=lnc --stranded=TRUE --outprefix=COUNT_FILE.lncRNA."
  o5 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste(cmd5.1,cmd5.2,sep="\n"),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)),
                                                                        pattern="COUNT_FILE",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/htseqCounts/",filename)))))
  write(o5, file="~/sandbox/htseqCounts")
  scpFile(file.local="~/sandbox/htseqCounts", dir.remote="~/bin/")
  
  
  
  
  
  
  
  cmd3 <- "java -jar -Xmx70g /share/pkg/picard/1.96/SamFormatConverter.jar INPUT=test.star.samAligned.out.sam OUTPUT=test.starPic.bam MAX_RECORDS_IN_RAM=5000000"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  write(o3, file="~/sandbox/starBamPicard")
  scpFile(file.local="~/sandbox/starBamPicard", dir.remote="~/bin/")
  df.comb$picSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".starPic_bam"))
  
  
 "samtools sort -m 171798691840 test.star.sam test.star_sort"
  
  
  
  
  o <- paste0(getBedopsIntersect(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),tag="star"))
  
  write(o, file="~/sandbox/starToBed.sh")
  # cat ~/bin/starToBed.sh | xargs -I{}  perl ~/bin/runJob.pl -c 2 -m 20480 -W 600 -Q short -t "star2Bed" -i "{}"
  
  
  
  
  
  
  
  
  
  o <- paste0(getBedopsIntersect(file.path(rnaseqdir,"segemehl",df.comb$bare),tag="seg"))
  
  write(o, file="~/sandbox/segToBed.sh")
  #  cat ~/bin/segToBed.sh | xargs -I{}  perl ~/bin/runJob.pl -c 2 -m 20480 -W 600 -Q short -t "seg2Bed2" -i "{}"
  
  
 # ../gencodeV19.exonMerge.bed > ../gencodeV19.intronExDiff.bed
  
  o <- paste0(getBedopsIntersectStartWithBed(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),tag="star"))
  
  write(o, file="~/sandbox/starBedToElem.sh")
  
  o <- sortStarBed(file.path(rnaseqdir,"starSpikeIn",df.comb$bare))
  write(o, file="~/sandbox/starBedSort.sh")
  
  
  o <- getBedopsIntersect(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),tag="star")
  write(o, file="~/sandbox/getStarCts.sh")
  # cat ~/bin/getStarCts.sh | xargs -I{}  perl ~/bin/runJob.pl -c 4 -m 20480 -W 600 -Q short -t "starCts" -i "{}"
  
  
  notdone <- sapply(paste0(file.path(rnaseqdir,"starSpikeIn",df.comb$bare),".star_sort.bed"),function(x)hpc.file.exists(x))
  o <- sortStarBed(file.path(rnaseqdir,"starSpikeIn",df.comb$bare))[which(notdone == FALSE)]
  write(o, file="~/sandbox/starSortLongQ.sh")
  scpFile(file.local="~/sandbox/starSortLongQ.sh", dir.remote="~/bin/")
  # cat ~/bin/starSortLongQ.sh | xargs -I{}  perl ~/bin/runJob.pl -c 4 -m 20480 -W 6000 -Q long -t "longQstarBed" -i "{}"

  
  
}

# samtools view -bS test.seg.sam > test.seg.bam;;
#   bedtools bamtobed -i test.seg.bam > test.seg.bed


downloadCellDataLpa <- function(){
  
  
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  
  df.cytNuc.celltypes <- unique(subset(df, type=="fastq" & (localization == "nucleus" | localization == "cytosol"))[["cell"]])
  df.fastq <- subset(df, type == "fastq" & (localization == "cell") & (cell %in% df.cytNuc.celltypes))
  filenames <- c(df.cytNuc.restFastq$filename)
  remote.site <- cshl.rnaseq.dir
  rna.remote <- file.path("/project/umw_zhiping_weng/wespisea/","rna-seq/")
  web.file <- paste0(remote.site,filenames)
  remote.file <- paste0(rna.remote,filenames)
  
  len <- length(remote.file)
  split <- floor(len/2)
  
  
  md <- paste0(paste0("wget --continue ", web.file[1:split], " -O ", remote.file[1:split] ," &\n",
                      "wget --continue ", web.file[(split+1):len], " -O ", remote.file[(split +1):len] ),collapse= "; \n")
  write(md, file="~/sandbox/wgetCellPapPam.sh")
  scpFile(file.local="~/sandbox/wgetCellPapPam.sh", dir.remote="~/bin/")
  
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  
  
  o <- genStarAlignCmdSpikeIN.paramFile(rd1=file.path(rnaseqdir,df.comb$read1.filename),
                                        rd2=file.path(rnaseqdir,df.comb$read2.filename),
                                        outfile=file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star.sam")))
  write(o, file="~/sandbox/runStar.sh")
  scpFile(file.local="~/sandbox/runStar.sh", dir.remote="~/bin/")
  # cat ~/bin/runStar.sh | xargs -I{}  perl ~/bin/runJob.pl -c 16 -m 3072 -W 600 -Q short -t "runStar" -i "{}"
  df.comb$starAln <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star.samAligned.out.sam"))
  
  
  cmd1 <- "samtools view -bS test.star.samAligned.out.sam -o test.star.bam;;samtools sort -m 171798691840 test.star.bam test.star_sort"
  o1 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd1,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  #write(o1, file="~/sandbox/starBamSamtools.sh")
  #scpFile(file.local="~/sandbox/starBamSamtools.sh", dir.remote="~/bin/")
  df.comb$samtoolsSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam"))
  cmd3 <- "samtools index test.star_sort.bam"
  o3 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd3,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  write(o3, file="~/sandbox/stIdxStar")
  scpFile(file.local="~/sandbox/stIdxStar", dir.remote="~/bin/")
  df.comb$bamBai <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.bam.bai"))
  # cat ~/bin/stIdxStar | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 4096 -W 600 -Q short -t bamIndexTest -i "{}"
  
  
  
  cmd4 <- "samtools view -h test.star_sort.bam -o test.star_sort.sam"
  o4 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=cmd4,pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)))))
  #write(o4, file="~/sandbox/sortBam2Sam")
  #scpFile(file.local="~/sandbox/sortBam2Sam", dir.remote="~/bin/")
  df.comb$samSort <- file.path(rnaseqdir,"starSpikeIn",paste0(df.comb$bare,".star_sort.sam"))
  # catMissingLast.sh ~/bin/sortBam2Sam | xargs -I{} perl /home/aw30w/bin/runJob.pl -c 4 -m 8192 -W 600 -Q short -t sortBam2sam -i "{}"
  
  
  
  
  cmd6.1 <- paste("cat /project/umw_zhiping_weng/wespisea/flux-capicitor/v19annotParamsNIST14 >> test.flux.params && echo STDOUT_FILE test.flux.out >> test.flux.params &&",
                  "echo STDERR_FILE test.flux.err >> test.flux.params && echo STATS_FILE test.flux.stats >> test.flux.params && echo COVERAGE_FILE test.flux.coverage ",
                  ">> test.flux.params")
  cmd6.2 <- paste("/home/aw30w/bin/flux-capacitor-1.6.1/bin/flux-capacitor --threads 24 -p test.flux.params",
                  "-i  xxTESTINPUTxx.star_sort.bam -o test.flux.output")
  
  o6 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste0(cmd6.1,";;",cmd6.2),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14",filename)),
                                                                        pattern="xxTESTINPUTxx",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/",filename)  ))))
  #write(o6, file="~/sandbox/fluxCap")
  #scpFile(file.local="~/sandbox/fluxCap", dir.remote="~/bin/")
  # perl /home/aw30w/bin/runTask.pl -f ~/bin/fluxCap -c 24 -m 1024 -W 600 -Q short -t flux
  
  
  
  
  # /home/aw30w/bin/sjcount-master/sjcount -bam bam_file [-ssj junctions_output] [-ssc boundaries_output] [-log log_file]
  cmd7 <- paste0("/home/aw30w/bin/sjcount-master/sjcount3   -quiet  -bam  xxTESTINPUTxx.star_sort.bam -ssj test.ssj -ssc test.scc -log test.log" )
  o7 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste0(cmd7),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn/ssjCount",filename)),
                                                                        pattern="xxTESTINPUTxx",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/",filename)  ))))
  write(o7, file="~/sandbox/ssjcountRun")
  scpFile(file.local="~/sandbox/ssjcountRun", dir.remote="~/bin/")
  
  
  o.comb <- paste0(o1,";;",o3,";;",o4,";;",o6," &;;",o7)
  write(o.comb, file="~/sandbox/processStarAlign")
  scpFile(file.local="~/sandbox/processStarAlign", dir.remote="~/bin/")
  # perl /home/aw30w/bin/runTask.pl -f ~/bin/processStarAlign -c 26 -m 4096 -W 10:0 -Q short -t flux
  
  
  
  
  
  
  
  
  
  
  
  cmd5.1 <- "bedtools multicov -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.bed > COUNT_FILE.multicov.allTrans.uniqReads"
  cmd5.2 <- "bedtools multicov -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.singleTrans.bed > COUNT_FILE.multicov.singleTrans.uniqReads"
  cmd5.3 <- "bedtools multicov -D -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.bed > COUNT_FILE.multicov.allTrans.dupReads"
  cmd5.4 <- "bedtools multicov -D -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.intEx.lncPc.singleTrans.bed > COUNT_FILE.multicov.singleTrans.dupReads"
  cmd5.5 <- "bedtools multicov  -split -s  -bams test.star_sort.bam -bed /project/umw_zhiping_weng/wespisea/gtf/gencode.v19.grace.bed > COUNT_FILE.multicov.grace"
  
  o5 <- as.character(unlist(sapply(df.comb$bare, function(filename)gsub(x=gsub(x=paste(cmd5.1,cmd5.2,cmd5.3,cmd5.4,cmd5.5,sep="\n"),pattern="test", replacement=file.path(rnaseqdir,"starSpikeIn",filename)),
                                                                        pattern="COUNT_FILE",
                                                                        replacement=file.path(rnaseqdir,"starSpikeIn/multicovCounts/",filename)))))
  write(o5, file="~/sandbox/multicovCounts")
  scpFile(file.local="~/sandbox/multicovCounts", dir.remote="~/bin/")
  df.comb$multicovBare <- file.path(rnaseqdir,"starSpikeIn/multicovCounts/",df.comb$bare)
  
  
}


