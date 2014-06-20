#  cat gencode.v19.annotation.gtf | awk '{print $1,$2,$3,$4,$6,$7,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}'|sed 's/\"//g'|sed 's/;//g'

v19 <- "/home/wespisea/data/gencode.v19.annotation.gtf"
v19.lnc <- "/home/wespisea/data/gencode.v19.long_noncoding_RNAs.gtf"

v19.lncPc <- "/home/wespisea/data/gencode.v19.lncNoCodeMultiEx.pc.gtf"
v19.lncPc.sort <- "/home/wespisea/data/gencode.v19.lncNoCodeMultiEx.pc.sort.gtf"

extractGeneTrans <- function(line){
  paste0(strsplit(x=line ,split='\"')[[1]][c(2,4)],collapse=" ")
}

lnc.cpat <- "/home/wespisea/data/gencode.v19.lncRNA_transcripts.cpat"
pc.cpat <- "/home/wespisea/data/gencode.v19.pc_transcripts.cpat"

plotDir <- "/home/wespisea/work/research/researchProjects/coexpr/lncNET/plots/gencode/codingPotential/"

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
}
  
plotCPAT <- function(){
    comb <- getAllCPAT()
    comb.fr <- comb[c("transcript_id","gene_name","type","mRNA_size", "ORF_size", "Fickett_score", "Hexamer_score","coding_prob")]
    comb.fr$ORF_ratio <- comb.fr$ORF_size / comb.fr$mRNA_size
    
    comb.fr$log10_mRNA_size <- log10(comb.fr$mRNA_size)
    comb.fr$mRNA_size <- NULL
    
    comb.fr$log10_ORF_size <- log10(comb.fr$ORF_size)
    comb.fr$ORF_size <- NULL
    m <- melt(comb.fr,id.var=c("transcript_id","gene_name","type"))
    m$value <- as.numeric(m$value)
    
    m$biotype <- ifelse(m$type == "lnc", "lncRNA", "mRNA")
    ggplot(m, aes(x=value,fill=biotype))+geom_density(alpha=I(0.4))+
      theme_bw()+facet_wrap(~variable,scale="free")+
      ggtitle("GENCODE V19 -> analysis of coding potential\nlncRNA vs. mRNA distributions")+
      theme(legend.position="top")
    ggsave(paste(outdir,"CPAT-lncVmRNA-density.pdf"),height=7,width=8)
    
    ggplot(m, aes(x=type,y=value,fill=biotype))+geom_boxplot()+
      theme_bw()+facet_wrap(~variable,scale="free")+
      ggtitle("GENCODE V19 -> analysis of coding potential\nlncRNA vs. mRNA distributions")+
      theme(legend.position="top")
    ggsave(paste(outdir,"CPAT-lncVmRNA-boxplot.pdf"),height=7,width=8)
    
    ggplot(m, aes(x=value,fill=biotype,color=biotype))+stat_ecdf()+
      theme_bw()+facet_wrap(~variable,scale="free")+
      ggtitle("GENCODE V19 -> analysis of coding potential\nlncRNA vs. mRNA cumm. distributions")+
      theme(legend.position="top")
    ggsave(paste(outdir,"CPAT-lncVmRNA-cdf.pdf"),height=7,width=8)
    

    comb$kmeansGroup <- as.numeric(kmeans(scale(center=TRUE,scale=TRUE,x=as.matrix(comb[c("mRNA_size", "ORF_size","Fickett_score","Hexamer_score")])),2)$cluster)
    comb$groupType <- with(comb,paste(type,kmeansGroup))
    
    
    #pComb <- data.frame(PC1=p$x[,1],PC2=p$x[,2],fill=label,color=label)
    p <- prcomp(as.matrix(comb[c("mRNA_size", "ORF_size","Fickett_score","Hexamer_score")]),scale.=TRUE,center=TRUE)
    pComb <- data.frame(PC1=p$x[,1],PC2=p$x[,2],PC3=p$x[,3],type=comb$type,groupType=comb$groupType,gene_name=comb$gene_name)
    ggplot(pComb, aes(x=PC1,y=PC2, fill=type, color=groupType))+geom_point() +theme_bw()
    ggplot(pComb, aes(x=PC1, fill=type))+geom_density(alpha=I(0.3)) +theme_bw()
    
    ggplot(data.frame(cumStdev=c(0,cumsum(p$sdev))/sum(p$sdev),comp=0:4), aes(x=comp,y=cumStdev)) + 
      geom_line()+ylim(0,1)+
      theme_bw()+
      ggtitle("PCA on CPAT\nFeatures:\nmRNA_size, ORF_size, Fickett_score, Hexamer_score")
    ggsave(paste(outdir,"CPAT-cumVariance-components.png"),height=7,width=6)
    
    ggplot(comb, aes(x=coding_prob,fill=groupType))+geom_density(alpha=I(0.4))+theme_bw()+
      facet_grid(type~.,scale="free")
   
    ggplot(pComb, aes(x=PC1,y=PC2, fill=type, color=groupType,label=gene_name))+
      geom_point(alpha=I(0.9)) +theme_bw()+geom_rug()+
      xlim(-5,10)+ylim(-10,5)+
      scale_color_manual(values=c("red","green","blue","black")) + 
      facet_grid(type~.)+
      ggtitle("PCA components\nColor by biotype/kmeans(k=2) groups")
    ggsave(paste(outdir,"CPAT-c1c2-lncVmRNA-components.png"),height=7,width=6)
    
    ggplot(pComb, aes(x=PC1,y=PC3, fill=type, color=groupType,label=gene_name))+
      geom_point(alpha=I(0.9)) +theme_bw()+geom_rug()+
      scale_color_manual(values=c("red","green","blue","black")) + 
      facet_grid(type~.)+
      xlim(-5,10) +ylim(-5,3)
      ggtitle("PCA components\nColor by biotype/kmeans(k=2) groups")
    ggsave(paste(outdir,"CPAT-c1c3-lncVmRNA-components.png"),height=7,width=6)
    
    ggplot(pComb, aes(x=PC2,y=PC3, fill=type, color=groupType,label=gene_name))+
      geom_point(alpha=I(0.9)) +theme_bw()+geom_rug()+
      xlim(-10,5) +ylim(-5,3)+  
      scale_color_manual(values=c("red","green","blue","black")) + 
      facet_grid(type~.)+
      ggtitle("PCA components\nColor by biotype/kmeans(k=2) groups")
    ggsave(paste(outdir,"CPAT-c2c3-lncVmRNA-components.png"),height=7,width=6)
    
    
    
    
    cutoff <- seq(from=0,to=1,by=0.01)
    foundGenes <- sapply(cutoff,function(x)length(unique(sort(lncPcCPAT[which(lncPcCPAT$type == "lnc" & lncPcCPAT$coding_prob < x),"gene_name"]))))
    c.df <- data.frame(cutoff=cutoff,recoveredGenes=foundGenes)
    c.df$percentGeneRecovered <- 100* c.df$recoveredGenes/max(c.df$recoveredGenes)
    c.df$label <- " "
    c.df[1:10 * 10 - 9,"label"] <- paste0(as.character(format(c.df[1:10 * 10 - 9,"percentGeneRecovered"],digits=2)),"%")
    
    ggplot(c.df, aes(x=cutoff,y=recoveredGenes,label=label))+geom_line()+
      theme_bw()+ggtitle("CPAT coding prob. cutoff vs.\nnumber of remaining lncRNA genes")+
      geom_hline(slope=0,yintercept=c.df[which(c.df$cutoff == 0.25 ),"recoveredGenes" ],color="green")+
      geom_text(size=4)
    ggsave(paste(outdir,"CPAT-cutoffFoundGenes.pdf"),height=5,width=5)
    
    ggplot(c.df, aes(x=cutoff,y=percentGeneRecovered))+geom_line()+ylim(0,100)+
      theme_bw()+ggtitle("CPAT coding prob. cutoff vs.\nnumber of remaining lncRNA genes")+
      geom_hline(slope=0,yintercept=c.df[which(c.df$cutoff == 0.25 ),"percentGeneRecovered" ],color="green")+
      geom_text()
    ggsave(paste(outdir,"CPAT-cutoffFoundGenes-percent.pdf"),height=5,width=5)
    
    library(scatterplot3d)
    
    vecIn <- function(vec,range){
      (vec > range[1]) & (vec < range[2])
    }
    
   # attach(pComb)
    with(pComb[which(vecIn(pComb$PC1,c(-5,10)) & vecIn(pComb$PC2 , c(-10,5)) & vecIn(pComb$PC3,c(-5,3)) ),],
         scatterplot3d(PC1,PC3,PC2, 
                  type="h", main="3D Scatterplot",color=as.numeric(factor(type)), pch = 19))
   # detach(pComb)
    
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
  
  # start with all lncRNA in biotypes
  lncGenesTrans <- df[which(df$gene_type %in% c("processed_transcript","lincRNA","antisense") &
                              df$transcript_type %in% c("processed_transcript","lincRNA","antisense")),]
  
  transcriptBlackList <- c("snoRNA", "miRNA", "retained_intron", "sense_intronic","rRNA", "snRNA")
  # remove transcripts by biotypes we don't want
  lncGeneTransNotBlocked <- lncGenes[-which(lncGenes$transcript_type %in% transcriptBlackList),]
  lncPcCPAT <- getAllCPAT()
  cutoff <- 0.25
  #lncRNA transcripts below cutoff
  nonCodingLncTrans <- lncPcCPAT[which(lncPcCPAT$type == "lnc" & lncPcCPAT$coding_prob < cutoff),"transcript_id"]
  lncGeneTransClean <- lncGeneTransNotBlocked[which(lncGeneTransNotBlocked$transcript_id %in% nonCodingLncTrans),]
  
  exonsPerTransLnc <- as.data.frame(group_by(lncGeneTransClean,transcript_id) %.% 
                               filter(feature=="exon") %.% 
                               summarise(exon_count = length(geneTrans_type),
                                         trans_type = paste0(transcript_type,collapse=" "),
                                         gene_type = gene_type[1]))
  
  # get a list of all the single exon transcripts
  ggplot(exonsPerTransLnc,aes(x=exon_count,width=0.1,fill=gene_type)) + 
    geom_bar(stat="bin",binwidth=1) + 
    xlim(0,15) + 
    ggtitle("Distribution of exons per transcript\nFill by gene_type")
  
  
  singleExonTranscripts <- exonsPerTransLnc[which(exonsPerTransLnc$exon_count < 2),"transcript_id"]
  
  
  pcGenes <- df[which(df$gene_type == "protein_coding"),]
  pcGenes$gtid <- with(pcGenes,paste0(gene_id," ",transcript_id))
  
  # remove single exon transcripts from lncs
  lncGeneTransCleanMultiExon <- lncGeneTransClean[which(!lncGeneTransClean$transcript_id %in% singleExonTranscripts),]
  
  # get all genes with at least one transcript
  cme.hasTrans <- unique(lncGeneTransCleanMultiExon[which(lncGeneTransCleanMultiExon$feature == "transcript"),"gene_id"])
  #finally, remove genes without at least one transcript(dangling genes)
  lncGeneTransCleanMultiExon.hasTrans <- lncGeneTransCleanMultiExon[which(lncGeneTransCleanMultiExon$gene_id %in% cme.hasTrans),]
  
  lncGeneTransCleanMultiExon.hasTrans$gtid <- with(lncGeneTransCleanMultiExon.hasTrans,paste0(gene_id," ",transcript_id))
  all.gtid <- unique(c(lncGeneTransCleanMultiExon.hasTrans$gtid,pcGenes$gtid))
  
  v19.out <- v19.txt[which(as.character(v19.txt.geneTrans) %in% all.gtid)]
  t <- tempfile()
  writeLines(v19.out,v19.lncPc)
  cmd <- paste("sort -T /project/umw_zhiping_weng/wespisea/tmp -k1,1V -k2,2g",v19.lncPc," > ",v19.lncPc.sort)
  system(cmd)

}


