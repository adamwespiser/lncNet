 library(IRanges) 
 library(GenomicRanges) 
 library(Rsamtools) 
 #source("http://bioconductor.org/biocLite.R")
 ##biocLite("IRanges")
 #biocLite("GenomicRanges")
 #biocLite("Rsamtools")
 
 #http://www.nbic.nl/uploads/media/Practical_exon_and_gene_quantification_in_R.pdf
 
fname <- "./data//gencode.v19.long_noncoding_RNAs.gtf"
regions <- read.table(fname, sep="\t", header=FALSE)
k <- which ( regions[,3] == 'exon')
exons <- regions[k,]
 
k <- which( levels(exons[,1])[exons[,1] ] == "chr1" ) 
exons.chr1 <- exons[k, ] 
 
 lst.exons.chr1 <- lapply( unique(exons.chr1[,1]), function(x, D){ 
   i <- which(D[,1] == x); 
   IRanges( start=D[i,4], end=D[i,5], names=D[i,9]) ; 
 }, D=exons.chr1 )
 
 names(lst.exons.chr1) <- unique(exons.chr1[,1])
 
 rl.exons.chr1 <- RangesList() 
 
 for( x in names(lst.exons.chr1) ) { 
   rl.exons.chr1[[x]] <-  lst.exons.chr1[[x]] } 
 
 
 
 