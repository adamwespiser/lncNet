

spikeInFile <<- getFullPath("./data/SpikeIn.tab")



getSpikeInDf <- function(){
  df <- read.csv(file = spikeInFile, header=TRUE, stringsAsFactors=FALSE, sep=",")
  df$erccString <- as.character(df$ERCC)
  df$gene_id <- paste0("ERCC-",str_replace_all(sprintf("%5i", df$ERCC), " ", "0"),"_gene")
  df
}

