

spikeInFile <<- getFullPath("./data/SpikeIn.tab")



getSpikeInDf <- function(){
  df <- read.csv(file = spikeInFile, header=TRUE, stringsAsFactors=FALSE, sep=",")
  df$erccString <- as.character(df$ERCC)
  df$gene_id <- paste0("ERCC-",str_replace_all(sprintf("%5i", df$ERCC), " ", "0"),"_gene")
  df
}


getEncodeSpikeInConc <- function(){
  spike.data <- getSpikeInDf()
f <- getFullPath("data/wgEncodeChslLongSpikeinsConcentrations.tsv")
ds <- read.csv(file=f, sep="\t", comment.char="#", stringsAsFactors=FALSE,header=FALSE)
colnames(ds)<-c("gene_id_short", "concentration_ENCODE")
ds$ERCC_number <- as.numeric(sapply(ds$gene_id_short, function(x)gsub(x=x,pattern="ERCC-", replacement="")))
ds$gene_id <- paste0(ds$gene_id_short, "_gene")  
spike.m <- merge(ds,spike.data,by="gene_id")
g <- ggplot(spike.m, aes(concentration_ENCODE, Pool14nmol.ul))+geom_point()
spike.m
}
  