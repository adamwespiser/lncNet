




convertTransToTransFull <- function(transFile,geneFile){
  tf <- tempfile()
  system( paste("cat",transFile," | sed 's/[;\"]//g' |awk -F' ' '{print $10,$12,$14}' > ",tf))
  trans.df <- read.csv(file=tf, sep=" ", stringsAsFactors=FALSE,header=FALSE)
  file.remove(tf)
  colnames(trans.df) <- c("gene_id", "transcript_id", "RPKM")
  trans.df$RPKM = as.numeric(trans.df$RPKM)
  exportAsTable(df=trans.df  ,file=geneFile)
}


processCellsMaxTransExpr <- function(){
  annot.df <- getRpkmFromBamDataForOneCell()
  annot.df <- annot.df[which(annot.df$rnaExtract == "longPolyA"),]
  annot.df <- annot.df[-which(annot.df$cell == "H1-hESC"),]
  genes <- annot.df$rfbGene
  trans <- annot.df$rpkmFromBamFile
  annot.df$rep <- ifelse(annot.df$replicate >2,annot.df$replicate -2,annot.df$replicate )
  annot.df$transFullRPKM <- gsub(x=annot.df$rfbGene,pattern="genes",replacement="transFullRPKM")
  transFullRPKM <- annot.df$transFullRPKM 
  sapply(seq_along(transFullRPKM), function(x)convertTransToTransFull(transFile=trans[x],geneFile=transFullRPKM[x]))
  df.together <- data.frame()
  for ( cell in unique(annot.df$cell)){
    a.cell <- annot.df[which(annot.df$cell == cell),]
    a.cell.cyt1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "cytosol"),"transFullRPKM"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt1$loc <- "cyt"
    a.cell.cyt1$rep <- 1
    a.cell.nuc1 <- read.csv(file=a.cell[which(a.cell$rep == 1 & a.cell$localization == "nucleus"),"transFullRPKM"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc1$loc <- "nuc"
    a.cell.nuc1$rep <- 1
    a.cell.cyt2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "cytosol"),"transFullRPKM"],sep="\t",stringsAsFactors=FALSE)
    a.cell.cyt2$loc <- "cyt"
    a.cell.cyt2$rep <- 2
    a.cell.nuc2 <- read.csv(file=a.cell[which(a.cell$rep == 2 & a.cell$localization == "nucleus"),"transFullRPKM"],sep="\t",stringsAsFactors=FALSE)
    a.cell.nuc2$loc <- "nuc"
    a.cell.nuc2$rep <- 2
    
    comb <- rbind(a.cell.cyt1,a.cell.cyt2,a.cell.nuc1,a.cell.nuc2)
    transExpr <- as.data.frame(group_by(comb,gene_id,transcript_id) %.% summarise(sum(RPKM)))
    colnames(transExpr) <- c("gene_id","transcript_id","RPKMsum")
    transExpr$RPKM_tieBreaker <- transExpr$RPKMsum + runif(seq_along(transExpr$RPKMsum))/(10^9)
    gene.df <- as.data.frame(group_by(transExpr, gene_id) %.% filter(RPKM_tieBreaker == max(RPKM_tieBreaker))) 
    cellTranscripts <- gene.df$transcript_id
    cTrans <- comb[which(comb$transcript_id %in% cellTranscripts),]
    cTrans$cell <- cell                       
    if (i == 1){
      df.together <- cTrans
      
    } else{
      df.together <- rbind(df.together,cTrans)
    }
  }
  exportAsTable(file=getFullPath("/data/rpkmFromBam-TopTransCellType.tab"),file=df.together)
}




getDataExprBoth <- function(){
  df.rfb.ratio <- readInTable(file=getFullPath("/data/rpkmFromBamCapData-lpa-cytFrac.tab"))
  df.rfb.ratio$cytFracPseudo <- with(df.lpa.ratio, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  df.rfb.ratio$cytFrac <- with(df.rfb.ratio, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  df.rfb.ratio.pos <- df.rfb.ratio[which(df.rfb.ratio$value.ave.cyt != 0 & df.rfb.ratio$value.ave.nuc != 0),]
  
  
  df.rfb.ratio.rpkm <- df.rfb.ratio.pos[which(df.rfb.ratio.pos$variable =="RPKM"),]
  
  
  df.flux.ratio <- readInTable(file=getFullPath("/data/fluxCapData-lpa-proc.tab"))
  df.flux.ratio$cytFracPseudo <- with(df.flux.ratio, (value.rep1.pseudo.cyt+value.rep2.pseudo.cyt)/(value.rep1.pseudo.cyt + value.rep2.pseudo.cyt + value.rep1.pseudo.nuc + value.rep2.pseudo.nuc))
  df.flux.ratio$cytFrac <- with(df.flux.ratio, (value.ave.cyt)/(value.ave.cyt + value.ave.nuc))
  
  df.flux.ratio.pos <- df.flux.ratio[which(df.lpa.ratio$value.ave.cyt != 0 & df.flux.ratio$value.ave.nuc != 0),]
  
  
  df.flux.ratio.rpkm <- df.flux.ratio.pos[which(df.flux.ratio.pos$variable =="RPKM"),]
  
}