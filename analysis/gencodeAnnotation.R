#  cat gencode.v19.annotation.gtf | awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\"//g'|sed 's/;//g'

v19 <- "/home/wespisea/data/gencode.v19.annotation.gtf"
v19.lnc <- "/home/wespisea/data/gencode.v19.long_noncoding_RNAs.gtf"
extractGeneTrans <- function(line){
  paste0(strsplit(x=line ,split='\"')[[1]][c(2,4)],collapse=" ")
}

lnc.cpat <- "/home/wespisea/data/gencode.v19.lncRNA_transcripts.cpat"
pc.cpat <- "/home/wespisea/data/gencode.v19.pc_transcripts.cpat"

readInCPAT <- function(biotype="lnc"){
  stopifnot(biotype %in% c("lnc","pc"))
  f <- ifelse(biotype == "lnc", lnc.cpat,pc.cpat)
  
  d <- read.csv(file=f,stringsAsFactors=FALSE,sep="\t")
  d.rows <- rownames(d)
  d$transcript_id <- sapply(d.rows, function(x)as.character(strsplit(x=x ,split='\\|')[[1]][1]))
  d$gene_id <- sapply(d.rows, function(x)as.character(strsplit(x=x ,split='\\|')[[1]][2]))
  d$gene_name <- sapply(d.rows, function(x)as.character(strsplit(x=x ,split='\\|')[[1]][6]))
  d$length <- sapply(d.rows, function(x)as.character(strsplit(x=x ,split='\\|')[[1]][7]))
  
  d$type <- biotype
  d
}

getAllCPAT <- function(plot=FALSE){
  lnc <- readInCPAT("lnc")
  pc <- readInCPAT("pc")
  pc$type = "mRNA"
  comb <- rbind(lnc,pc)
  if(plot){
    
    
    comb.fr <- comb[c("transcript_id","gene_name","type","mRNA_size", "ORF_size", "Fickett_score", "Hexamer_score","coding_prob")]
    comb.fr$ORF_ratio <- comb.fr$ORF_size / comb.fr$mRNA_size
    
    comb.fr$log10_mRNA_size <- log10(comb.fr$mRNA_size)
    comb.fr$mRNA_size <- NULL
    
    comb.fr$log10_ORF_size <- log10(comb.fr$ORF_size)
    comb.fr$ORF_size <- NULL
    m <- melt(comb.fr,id.var=c("transcript_id","gene_name","type"))
    m$value <- as.numeric(m$value)
    ggplot(m, aes(x=value,fill=type))+geom_density(alpha=I(0.4))+
      theme_bw()+facet_wrap(~variable,scale="free")
    
    
    ggplot(comb, aes(x=coding_prob,fill=type))+geom_density(alpha=I(0.4))+
      theme_bw()+ggtitle("CPAT coding probability\nmRNA & lncRNA")
    
    
    ggplot(comb, aes(x=coding_prob,fill=type))+geom_bar(binwidth=0.05)+
      theme_bw()+ggtitle("CPAT coding probability\nmRNA & lncRNA") +
      facet_grid(type~.)
    ggplot(comb, aes(x=coding_prob))+stat_ecdf()+
      theme_bw()+ggtitle("CPAT coding probability\nmRNA & lncRNA\ncummulative prob.") +
      facet_grid(type~.)
    #ORF_size
    ggplot(comb, aes(x=log10(ORF_size)))+geom_freqpoly()+
      theme_bw()+ggtitle("CPAT ORF_size\nmRNA & lncRNA\nlog10 distribution") +
      facet_grid(type~.)
    
    ggplot(comb, aes(x=log10(ORF_size),y=log10(mRNA_size))) + 
      geom_rug(alpha=I(0.3))+ geom_point(size=1)+ geom_density2d()+
      xlab("log10(ORF length)")+ylab("log10(transcript length)")+
      theme_bw()+ ggtitle("CPAT ORF_size,mRNA_size\nmRNA & lncRNA\nlog10 distribution") +
      facet_grid(type~.)
    
    
    comb$kmeansGroup <- as.numeric(kmeans(scale(center=TRUE,scale=TRUE,x=as.matrix(comb[c("mRNA_size", "ORF_size","Fickett_score","Hexamer_score")])),2)$cluster)
    comb$groupType <- with(comb,paste(type,kmeansGroup))
    
    
    #pComb <- data.frame(PC1=p$x[,1],PC2=p$x[,2],fill=label,color=label)
    p <- prcomp(as.matrix(comb[c("mRNA_size", "ORF_size","Fickett_score","Hexamer_score")]),scale.=TRUE,center=TRUE)
    pComb <- data.frame(PC1=p$x[,1],PC2=p$x[,2],type=comb$type,groupType=comb$groupType,gene_name=comb$gene_name)
    ggplot(pComb, aes(x=PC1,y=PC2, fill=type, color=groupType))+geom_point() +theme_bw()
    ggplot(pComb, aes(x=PC1, fill=type))+geom_density(alpha=I(0.3)) +theme_bw()
    
    
     
    ggplot(comb, aes(x=coding_prob,fill=groupType))+geom_density(alpha=I(0.4))+theme_bw()+
      facet_grid(type~.,scale="free")
    ggplot(pComb, aes(x=PC1,y=PC2, fill=type, color=groupType,label=gene_name))+
      geom_point(alpha=I(0.9)) +theme_bw()+
      xlim(-5,10)+ylim(-10,5)+geom_text(size=3)+
      scale_color_manual(values=c("red","green","blue","black")) + 
      facet_grid(type~.)
    
    
    
  }
}


readInGtfGencode <- function(file = v19){
  gtfCols <- c("chr", "annot", "feature", "start", "uk1", "uk2","uk3", "gene_id", "transcript_id", 
               "gene_type", "gene_status", "gene_name", "transcript_type","transcript_status","transcript_name",
               "exon_number","exon_id", "havana_gene","havana_transcript")
  
  p1 <- pipe(paste("cat",v19,"| awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\\\"//g'|sed 's/;//g'" ))
  v19.lines <- readLines(v19)
  v19.txt <- v19.lines[-grep(x=v19.lines,pattern="#")]
  v19.txt.geneTrans <- sapply(v19.txt,extractGeneTrans)
  df <- read.csv(file=p1,sep=" ",stringsAsFactors=FALSE,comment.char="#",header=FALSE)
  colnames(df) <- gtfCols
  
  p2.lnc <- pipe(paste("cat",v19.lnc,"| awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\\\"//g'|sed 's/;//g'" ))
  df.lnc <- read.csv(file=p2.lnc,sep=" ",stringsAsFactors=FALSE,comment.char="#",header=FALSE)
  colnames(df.lnc) <- gtfCols
  
  df.lnc$geneTrans_status <- with(df.lnc,paste0(gene_status,"--",transcript_status))
  df.lnc$geneTrans_type <- with(df.lnc,paste0(gene_type,"--",transcript_type))
  

  
  
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
  lncTrans <- df[which(df$transcript_type %in% c("processed_transcript","lincRNA","antisense")),]
  exonsPerT <- as.data.frame(group_by(,gene_id) %.% 
                                   filter(feature=="transcipt") %.% 
                                   summarise(type_count = length(unique(geneTrans_type))))
  
  lncGenesTrans <- df[which(df$gene_type %in% c("processed_transcript","lincRNA","antisense") &
                              df$transcript_type %in% c("processed_transcript","lincRNA","antisense")),]
  
  transcriptBlackList <- c("snoRNA", "miRNA", "retained_intron", "sense_intronic","rRNA", "snRNA")
  lncGeneTransClean <- lncGenesTrans[-which(lncGenesTrans$transcript_type %in% transcriptBlackList),]
  exonsPerTransLnc <- as.data.frame(group_by(lncGeneTransClean,transcript_id) %.% 
                               filter(feature=="exon") %.% 
                               summarise(exon_count = length(geneTrans_type),
                                         trans_type = paste0(transcript_type,collapse=" "),
                                         gene_type = gene_type[1]))
  
  ggplot(exonsPerTransLnc,aes(x=exon_count,width=0.1,fill=gene_type)) + 
    geom_bar(stat="bin",binwidth=1) + 
    xlim(0,15) + 
    ggtitle("Distribution of exons per transcript\nFill by gene_type")
  
  pcGenes <- df[which(df$gene_type == "protein_coding"),]
  pcGenes$gtid <- with(pcGenes,paste0(gene_id," ",transcript_id))
  
  singleExonTranscripts <- exonsPerTransLnc[which(exonsPerTransLnc$exon_count == 1),"transcript_id"]
  lncGeneTransCleanMultiExon <- lncGeneTransClean[which(!lncGeneTransClean$transcript_id %in% singleExonTranscripts),]
  cme.hasTrans <- lncGeneTransCleanMultiExon[which(lncGeneTransCleanMultiExon$feature == "transcript"),"gene_id"]
  lncGeneTransCleanMultiExon.hasTrans <- lncGeneTransCleanMultiExon[which(lncGeneTransCleanMultiExon$gene_id %in% cme.hasTrans),]
  lncGeneTransCleanMultiExon.hasTrans$gtid <- with(lncGeneTransCleanMultiExon.hasTrans,paste0(gene_id," ",transcript_id))
  all.gtid <- unique(c(lncGeneTransCleanMultiExon.hasTrans$gtid,pcGenes$gtid))
  
  v19.out <- v19.txt[which(as.character(v19.txt.geneTrans) %in% all.gtid)]

}


