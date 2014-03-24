homeFolder <- path.expand("~")
getDataFile <- function(subpath){file.path(homeFolder, "data", subpath)}
projectDir <- normalizePath(file.path(script.dir(), ".."))
getFullPath <- function(subpath){ file.path(projectDir, subpath) }

source(getFullPath("analysis/gencodeInput.R"))


local.datadir <<- "/home/wespisea/data/"

gene.expr <<- getDataFile("gencodev10_long_genes_with_quantif_and_npIDR_83exp.txt")
gene.expr.tab <<- getDataFile("gencodev10_long_genes_with_quantif_and_npIDR_83expLong.tab")
gene.expr.annot <<- getDataFile("83exp.txt")
lnc.Annot <<- "/home/wespisea/work/research/researchProjects/encode/encode-manager/data/lncList_ensemblIdsFoundLncRna2.list"

getLncFunc <- function(f=lnc.Annot){
  read.csv(file=f, stringsAsFactors=FALSE, sep ="\t")
}


getAnnotFile <- function(f=gene.expr.annot){
  df <- read.csv(f,sep=".",header=FALSE,stringsAsFactors=FALSE)
  df$V1 <- sapply(df$V1, function(x)gsub(x,pattern=" rnafrac", replacement=""))
  df <- df[c(1,2,4,6)]
  colnames(df) <- c("expr", "pulldown", "cell", "localization")
  df$annot <- paste0(df$cell,".", df$localization, ".", df$pulldown, sep="")
  df
}
df <- getAnnotFile()

mapExprTNumeric <- function(exprStrVec,index,splitChar){
  sapply(levels(exprStrVec)[exprStrVec],  function(label)as.numeric(unlist(strsplit(label,splitChar)[[1]])[index]))
}
mapExprTChar <- function(exprStrVec,index,splitChar){
  sapply(levels(exprStrVec)[exprStrVec],  function(label)unlist(strsplit(label,splitChar)[[1]])[index])
}
mapExprTChar_char <- function(exprStrVec,index,splitChar){
  sapply(exprStrVec,  function(label)unlist(strsplit(label,splitChar)[[1]])[index])
}
processExprFile <- function(f=gene.expr){
  annot <- getAnnotFile()
  annot$colTitle <- sapply(annot$expr, function(x)gsub(x,pattern="-",replacement="\\."))
  df <- read.csv(file=gene.expr, sep=" ",stringsAsFactors=FALSE)
  annot$colTitle %in% names(df)[-1]
  newCols <- as.character(sapply(names(df)[-1], function(x)annot[which(annot$colTitle == x),"annot"]))
  colnames(df) <- c("gene_id",newCols)
  melt.df <- melt(df, id.var = "gene_id")
  splitExpr <- strsplit(levels(melt.df$value)[melt.df$value],"\\:")
  melt.df$RPKM1 <- sapply(splitExpr,function(x)x[1])
  melt.df$RPKM2 <- sapply(splitExpr,function(x)x[2])
  melt.df$IDR <- sapply(splitExpr,function(x)x[3])
  splitAnnot <- strsplit(levels(melt.df$variable)[melt.df$variable],"\\.")
  melt.df$cell <- sapply(splitAnnot,function(x)x[1])
  melt.df$localization <- sapply(splitAnnot,function(x)x[2])
  melt.df$pulldown <- sapply(splitAnnot,function(x)x[3])
  splitAnnot <- strsplit(melt.df$gene, "\\.")
  melt.df$gene_id_short <- sapply(splitAnnot,function(x)x[1])
  
  lnc.genes <- getGencodeAnnot("lnc", "v10")[["gene_id"]]
  pc.genes <- getGencodeAnnot("pc", "v10")[["gene_id"]]
  melt.df$biotype <- "other"
  melt.df[which(melt.df$gene_id %in% pc.genes),"biotype"] <- "pc"
  melt.df[which(melt.df$gene_id %in% lnc.genes),"biotype"] <- "lnc"
  
  lncFunc <- getLncFunc()
  melt.df$funcLnc <- 0
  melt.df[which(melt.df$gene_id_short %in% lncFunc$ensembl_gene_id),"funcLnc"] <- 1
  
  melt.df$RPKMsum = melt.df$RPKM1 + melt.df$RPKM2
  
  
  exportAsTable(df=melt.df,file=gene.expr.tab )
}

getExprFile <- function(f=gene.expr.tab){
  if(!file.exists(f)){
    processExprFile()
  }
  d <- read.csv(file=f, stringsAsFactors = FALSE, sep = "\t", 
                colClasses =c(rep("character",3), rep("numeric",3),rep("character",5),"numeric") )
  df$RPKMsum = df$RPKM1 + df$RPKM2
  
}

getCytNuc <- function(){
  annot <- getAnnotFile()
  annot.lpa <- subset(annot, pulldown == "longPolyA")
  annot.lnpa <- subset(annot, pulldown == "longNonPolyA")
  annot.lpa.cell <- subset(annot.lpa, localization == "cell") 
  annot.lpa.cyt <- subset(annot.lpa, localization == "cytosol")
  annot.lpa.nuc <- subset(annot.lpa, localization == "nucleus")
  annot.lnpa.cell <- subset(annot.lnpa, localization == "cell") 
  annot.lnpa.cyt <- subset(annot.lnpa, localization == "cytosol")
  annot.lnpa.nuc <- subset(annot.lnpa, localization == "nucleus")
  
  
}

