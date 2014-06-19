#  cat gencode.v19.annotation.gtf | awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\"//g'|sed 's/;//g'

v19 <- "/home/wespisea/data/gencode.v19.annotation.gtf"

readInGtfGencode <- function(file = v19){
  p1 <- pipe(paste("cat",v19,"| awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\\\"//g'|sed 's/;//g'" ))
  df <- read.csv(file=p1,sep=" ",stringsAsFactors=FALSE,comment.char="#",header=FALSE)
  colnames(df) <- c("chr", "annot", "feature", "start", "uk1", "uk2","uk3", "gene_id", "transcript_id", 
                    "gene_type", "gene_status", "gene_name", "transcript_type","transcript_status","transcript_name",
                    "exon_number","exon_id", "havana_gene","havana_transcript")
  df$geneTrans_status <- with(df,paste0(gene_status,"--",transcript_status))
  df$geneTrans_type <- with(df,paste0(gene_type,"--",transcript_type))
  typeKeep <- c("antisense--antisense","antisense--retained_intron",'lincRNA--lincRNA',
                "lincRNA--processed_transcript","processed_transcript--antisense",
                "processed_transcript--lincRNA","processed_transcript--processed_transcript",
                "processed_transcript--retained_intron","protein_coding--protein_coding")
  
  
  
  exonsPerGenes <- as.data.frame(group_by(df,gene_id) %.% 
                                   filter(feature=="exon") %.% 
                                   summarise(exonPerGene=length(uk1),
                                             gene_type = gene_type[1],
                                              type_count = length(unique(geneTrans_type))))
  
  lncGenes <- df[which(df$gene_type %in% c("processed_transcript","lincRNA","antisense")),]
  exonsPerT <- as.data.frame(group_by(,gene_id) %.% 
                                   filter(feature=="transcipt") %.% 
                                   summarise(type_count = length(unique(geneTrans_type))))
  
  
  
  
}