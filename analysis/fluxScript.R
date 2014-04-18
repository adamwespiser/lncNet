library(ggplot2)
library(dplyr)
library(plyr)
library(VennDiagram)

rnaseqdir <<- "/project/umw_zhiping_weng/wespisea/rna-seq/"

ghpc <<- "aw30w@ghpcc06.umassrc.org"
lnc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.long_noncoding_RNAs.geneList"
pc.v19.list <<- "/home/wespisea/data/gtf//gencode.v19.annotation.pc.geneList"



filterFluxCapictor <- function(fileIn,fileOut){
  paste("cat ",fileIn, "| sed 's/[;\\\"]//g' | awk -F ' ' '{if ($3 == \"transcript\"){print $3,$14,$16,$18,$20} if($3 == \"junction\"){len=$5-$4;print $3,$10,$14,len,\"0.0\";} else if ($3 == \"intron\"){len=$5-$4;print $3,$10,$14,len,\"0.0\"}}' >",
        fileOut)
}


thisTheme <<- theme_bw() +
  theme(text = element_text(size=12)) + 
  theme(panel.grid.major.x = element_line(colour = "grey"))+
  theme(panel.grid.minor.x = element_line(colour = "grey")) +
  theme(panel.grid.major.y = element_blank())+
  theme(panel.grid.minor.y = element_blank()) +
  theme(legend.position = "top")






readInFluxGtfParsed <- function(file="/home/wespisea/data/flux/wgEncodeCshlLongRnaSeqK562NucleusPapFastqRep2.flux.output.proc"){
  df <- read.csv(file=file,stringsAsFactors=FALSE, header=FALSE,sep=" ")
  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  colnames(df) <- c("region", "gene_id", "reads", "length", "RPKM_byFluxC")
  
  df$lengthKb <- df$length/1000
  df.intron <- df[which(df$region == "intron"), ]
  df.junction <- df[which(df$region == "junction"), ]
  df.transcript <- df[which(df$region == "transcript"), ]
  df.spike <- df.transcript[grep(pattern="ERCC",df.transcript$gene_id),]
  spikeInReads <- sum( df.spike$reads )
  spikeInReadsPerKb <- mean(df.spike$reads / df.spike$lengthKb) 
  
  totalReads <- sum(df.transcript$reads,na.rm=TRUE)
  millionsOfReads <- totalReads/(10^6)
  spikeInReadsPerKbMill <- spikeInReadsPerKb/millionsOfReads
  
  millReads <- max(with(df.transcript,1/( RPKM_byFluxC/(reads/(lengthKb)))),na.rm=TRUE)
  df.transcript$readsPerKb = with(df.transcript, (reads/lengthKb))
  df.transcript$RPKM = with(df.transcript, readsPerKb/millionsOfReads)
  df.transcript$RPKM_spikeIn = with(df.transcript, RPKM/spikeInReadsPerKbMill)
  df.transcript$readsPerLength = with(df.transcript, (reads/length))
  
  
  
  df.transcript$one <- 1
  df.gene.trans <- df.transcript %.% 
    group_by("gene_id") %.% 
    summarize(transTotalRPKM = sum(RPKM),
              transcriptTotalReadsPerKb = sum(readsPerKb),
              transcriptTotalRPKM_spikeIn = sum(RPKM_spikeIn),
              transcriptsPerGene = sum(one),
              transcriptTotalReads = sum(reads),
              transcriptMaxLength = max(length),
              transcriptAveLength = max(length)/sum(one),
              transcriptReadsPerLength = sum(readsPerLength))
  df.gene.trans$transcriptCount <- sum(ifelse(df.gene.trans$transcriptTotalReads > 0 & !is.na(df.gene.trans$transcriptTotalReads),1,0))
  df.junction$one <- 1
  df.junc.trans <- df.junction %.% 
    group_by("gene_id") %.% 
    summarize(junctionTotalReads = sum(reads), 
              junctionCount = sum(one) )
  
  
  df.intron$one <- 1
  df.intron$intronRPKM = with(df.intron, (reads/length)*millReads)
  df.intron$intronReadsPerLen = with(df.intron, (reads/length))
  df.intron.trans <- df.intron %.% 
    group_by("gene_id") %.% 
    summarize(intronTotalReads = sum(reads), 
              intronTotalLength = sum(length),
              intronCount = sum(one),
              #   intronTotalRPKM = sum(intronRPKM),
              intronReadsPerLength = sum(intronReadsPerLen))
  df.intron.trans$intronCount <- sum(ifelse(df.intron.trans$intronTotalReads > 0 & !is.na(df.intron.trans$intronTotalReads),1,0))
  df.intron.trans$intronTotalRPKM = (df.intron.trans$intronTotalReads/ ((df.intron.trans$intronTotalLength/1000) * millionsOfReads))
  df.intron.trans$intronTotalRPKM_spikeIn = df.intron.trans$intronTotalRPKM / millionsOfReads
  
  df.comb <- merge(x=merge(as.data.frame(df.gene.trans), as.data.frame(df.junc.trans),by="gene_id", all.x=TRUE ),as.data.frame(df.intron.trans), by="gene_id", all.x=TRUE)
  df.comb$gene_type <- "other"
  df.comb[which(df.comb$gene_id %in% unique(pc)),"gene_type"] <- "mRNA"
  df.comb[which(df.comb$gene_id %in% unique(lnc)),"gene_type"] <- "lncRNA"
  #df.comb$transTotalRPKMnorm      <- df.comb$transTotalRPKM  / mean(df.comb[grep(pattern="ERCC",df.comb$gene_id),"transTotalRPKM"])
  df.comb
}


plotIntronCount <- function(){
  df <- readInFluxGtfParsed()
  ggplot(df, aes(intronTotalReads/intronTotalLength, junctionTotalReads/junctionCount,color=gene_type))+geom_point()+
    theme_bw() + xlim(0,2) + ylim(0,1000) + scale_color_manual(values=c("red","green","yellow"))
  
  ggplot(df, aes(intronTotalReads/intronTotalLength, (junctionTotalReads/junctionCount),color=gene_type))+geom_point()+
    theme_bw() + xlim(0,2) + ylim(0,1000) + scale_color_manual(values=c("red","green","yellow"))
  
  ggplot(df, aes((intronTotalReads/intronCount), (junctionTotalReads/junctionCount),color=gene_type))+geom_point()+
    theme_bw()  + xlim(0,2000) + ylim(0,1000) + scale_color_manual(values=c("red","green","yellow"))
  
  ggplot(df, aes((intronTotalReads/intronTotalLength), (junctionTotalReads/junctionCount)/transTotalRPKM,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +xlim(0,2)
  
  ggplot(df, aes((intronTotalReads/intronTotalLength), transTotalRPKM,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) 
  
  ggplot(df, aes((intronTotalReads/intronTotalLength), transTotalRPKM,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) + xlim(0,5)
  
  ggplot(df, aes((intronTotalReads/intronTotalLength), transTotalRPKM,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) + xlim(0,2) + ylim(0,500)
  
  
  ggplot(df, aes((intronTotalLength)/intronCount, transcriptMaxLength,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) + 
    xlim(0,250000) + ylim(0,20000)
  
  ggplot(df, aes(transcriptAveLength,fill=gene_type))+geom_density(alpha=I(0.4)) + theme_bw() + xlim(0,2000)
  ggplot(df, aes(intronTotalLength,fill=gene_type))+geom_density(alpha=I(0.4)) + theme_bw() + xlim(0,200000)
  
  ggplot(df, aes(intronCount,junctionCount, color=gene_type))+geom_point() + theme_bw()
  
  ggplot(df[which(df$transTotalRPKM > 0),], aes(transTotalRPKM, fill =gene_type))+geom_density(alpha=I(0.4)) + theme_bw() +
    xlim(0,20)
  
  ggplot(df, aes(intronTotalRPKM, transTotalRPKM, color=gene_type)) + geom_point() +
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) + geom_abline(intercept=0,slope=1)
  
  ggplot(df, aes(intronTotalRPKM, transTotalRPKM, color=gene_type)) + geom_point() +
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) + geom_abline(intercept=0,slope=1) +
    xlim(0,200) + ylim(0,1500) +facet_wrap(~gene_type)
  
}


getFluxDataForOneCell <- function(){
  df <- read.csv(file=filesTxtTab, stringsAsFactors=FALSE, sep="\t")
  df.fastq <- subset(df,type=="fastq" & (localization == "nucleus" | localization == "cytosol"))
  read1 <- grep(df.fastq$filename,pattern="Rd1")
  read2 <- grep(df.fastq$filename,pattern="Rd2")
  
  
  df.comb <- data.frame(read1 = df.fastq[read1,], read2=df.fastq[read2,])
  df.comb$bare <- gsub(gsub(df.comb$read1.filename,pattern="Rd1",replacement=""),pattern=".fastq.gz",replacement="")
  df.comb$remote <- file.path(rnaseqdir,"starSpikeIn/flux-capacitorNIST14",paste0(df.comb$bare,".flux.output"))
  o1 <- paste0("scp aw30w@ghpcc06.umassrc.org:",df.comb$remote, " /home/wespisea/data/flux/")
  write(o1,file="~/sandbox/fetchFlux")
  o2 <- filterFluxCapictor(paste0("/home/wespisea/data/flux/",df.comb$bare,".flux.output"),
                           paste0("/home/wespisea/data/flux/",df.comb$bare,".flux.output.proc"))
  write(o2,file="~/sandbox/procFlux")
  
  df.comb$fluxFile <- paste0("/home/wespisea/data/flux/",df.comb$bare,".flux.output.proc")
  df.comb <- df.comb[c("read1.localization", "read1.cell", "read1.rnaExtract","read2.replicate" ,"fluxFile", "bare")]
  colnames(df.comb) <- c("localization", "cell", "rnaExtract","replicate" ,"fluxFile", "bare")
  df.comb
}



getTranscriptData <- function(celltype,rnaExtract){
  annot.df <- getFluxDataForOneCell()
  if (!missing(celltype)){
    annot.df <- annot.df[which(annot.df$cell == celltype),]
  }
  if(!missing(rnaExtract)){
    annot.df <- annot.df[which(annot.df$rnaExtract == rnaExtract),]
  }
  
  
  df.together <- data.frame()
  for(i in seq_along(annot.df$fluxFile)){
    df.local <- readInFluxGtfParsed(file=annot.df$fluxFile[i])
    df.local$cell <- annot.df$cell[i]
    df.local$localization <- annot.df$localization[i]
    df.local$rnaExtract <- annot.df$rnaExtract[i]
    df.local$replicate <- ifelse(annot.df$replicate[i] > 2, annot.df$replicate[i] -2, annot.df$replicate[i])
    if (i == 1){
      df.together <- df.local
      
    } else{
      df.together <- rbind(df.together,df.local)
    }
  }
  df.together
  
}



plotFluxData <- function(){
  df.together <- getTranscriptData()
  pc <- readLines(pc.v19.list)
  lnc <- readLines(lnc.v19.list)
  
  
  ggplot(df.together, aes(y=(intronTotalReads/intronTotalLength),x= transTotalRPKM,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +xlim(0,1000) + ylim(0,50)+
    ggtitle("K562 + GM12878 , all transcripts")
  gsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-intron-vs-exon-Reads.pdf"))
  
  ggplot(df.together[which(df.together$cell == "K562" & df.together$transTotalRPKMnorm <50 ),], aes(y=intronTotalRPKMnorm,x= transTotalRPKMnorm,color=gene_type))+geom_point()+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +
    facet_wrap(localization ~ rnaExtract,scale="free")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-intron-vs-exon-Reads-K562.pdf"))
  
  ggplot(df.together[which(df.together$cell == "K562"  ),], aes(x=log(intronTotalRPKMnorm),fill=gene_type))+geom_density(alpha=I(0.4))+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +#xlim(0,0.1)+
    facet_wrap(localization ~ rnaExtract) +
    ggtitle("intronRPKM / ave(spikeIn RPKM) : K562 ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-intronDensity-Normalized-K562.pdf"))
  
  ggplot(df.together[which(df.together$cell == "GM12878"  ),], aes(x=log(intronTotalRPKMnorm),fill=gene_type))+geom_density(alpha=I(0.4))+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +
    facet_wrap(localization ~ rnaExtract,scale="free") +
    ggtitle("intronRPKM / ave(spikeIn RPKM) : GM12878 ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-intronDensity-Normalized-GM12878.pdf"))
  
  ggplot(df.together[which(df.together$cell == "K562"  ),], aes(x=log(transTotalRPKMnorm),fill=gene_type))+geom_density(alpha=I(0.4))+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +
    facet_wrap(localization ~ rnaExtract,scale="free") +
    ggtitle("transRPKM / ave(spikeIn RPKM) : K562 ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-transDensity-Normalized-K562.pdf"))
  
  ggplot(df.together[which(df.together$cell == "GM12878"  ),], aes(x=log(transTotalRPKMnorm),fill=gene_type))+geom_density(alpha=I(0.4))+
    theme_bw() + scale_color_manual(values=c("red","green","yellow")) +#xlim(0,0.1)+
    facet_wrap(localization ~ rnaExtract ) +
    ggtitle("transRPKM / ave(spikeIn RPKM) : GM12878 ")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/flux-transDensity-Normalized-GM12878.pdf"))
  
  
  df.gm.lpa.cyt <- as.data.frame(group_by(df.together, "gene_id") %.%
                                   filter(cell == "GM12878" ) %.%
                                   filter(rnaExtract =="longPolyA") %.%
                                   filter(localization == "cytosol") %.%
                                   summarize(RPKMSum = sum(transTotalRPKMnorm)))
  
  df.gm.lpa.nuc <- as.data.frame(group_by(df.together, "gene_id") %.%
                                   filter(cell == "GM12878" ) %.%
                                   filter(rnaExtract =="longPolyA") %.%
                                   filter(localization == "nucleus") %.%
                                   summarize(RPKMSum = sum(transTotalRPKMnorm)))
  
  df.gm.lpa.ratio <- merge(df.gm.lpa.nuc,df.gm.lpa.cyt,by="gene_id",suffixes=c(".nuc",".cyt"))
  df.gm.lpa.ratio$RPKM.ratio <- with(df.gm.lpa.ratio, RPKMSum.cyt/RPKMSum.nuc)
  df.gm.lpa.ratio$gene_type <- "other"
  df.gm.lpa.ratio[which(df.gm.lpa.ratio$gene_id %in% lnc),"gene_type"] <- "lncRNA"
  df.gm.lpa.ratio[which(df.gm.lpa.ratio$gene_id %in% pc),"gene_type"] <- "mRNA"
  df.gm.lpa.ratio <- df.gm.lpa.ratio[which(df.gm.lpa.ratio$RPKMSum.nuc > 0 & df.gm.lpa.ratio$RPKMSum.cyt > 0 ),]
  
  
  
  
  ggplot(df.gm.lpa.ratio, aes(y=RPKMSum.cyt,x=RPKMSum.nuc)) + geom_point()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_abline(slope=1,intercept=0)+
    ggtitle("GM12878 normailzed RPKM, cyt, nuc")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratio-GM12878.pdf"))
  
  m <-  plot(df.gm.lpa.ratio, aes(y=log(RPKMSum.cyt),x=log(RPKMSum.nuc))) + geom_point() + 
    theme_bw() + geom_abline(slope=1,intercept=0)+
    ggtitle("GM12878 normailzed RPKM, cyt, nuc")+ xlim(0,4) + ylim(0,4)
  m + geom_density2d()
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratio-GM12878.pdf"))
  
  
  
  ggplot(df.gm.lpa.ratio, aes(y=RPKMSum.cyt,x=RPKMSum.nuc)) + geom_point()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_abline(slope=1,intercept=0) + 
    ggtitle("GM12878 normailzed RPKM, cyt, nuc") + xlim(0,10) + ylim(0,10)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratio-GM12878-zoom.pdf"))
  
  ggplot(df.gm.lpa.ratio, aes(log(RPKM.ratio))) + geom_density()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_vline(x=0)
  ggtitle("GM12878 normailzed RPKM, cyt / nuc") 
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratioDensity-GM12878-.pdf"))
  
  
  
  
  df.gm.lpa.cyt <- as.data.frame(group_by(df.together, "gene_id") %.%
                                   filter(cell == "K562" ) %.%
                                   filter(rnaExtract =="longPolyA") %.%
                                   filter(localization == "cytosol") %.%
                                   summarize(RPKMSum = sum(transTotalRPKMnorm)))
  
  df.gm.lpa.nuc <- as.data.frame(group_by(df.together, "gene_id") %.%
                                   filter(cell == "K562" ) %.%
                                   filter(rnaExtract =="longPolyA") %.%
                                   filter(localization == "nucleus") %.%
                                   summarize(RPKMSum = sum(transTotalRPKMnorm)))
  
  df.gm.lpa.ratio <- merge(df.gm.lpa.nuc,df.gm.lpa.cyt,by="gene_id",suffixes=c(".nuc",".cyt"))
  df.gm.lpa.ratio$RPKM.ratio <- with(df.gm.lpa.ratio, RPKMSum.cyt/RPKMSum.nuc)
  df.gm.lpa.ratio$gene_type <- "other"
  df.gm.lpa.ratio[which(df.gm.lpa.ratio$gene_id %in% lnc),"gene_type"] <- "lncRNA"
  df.gm.lpa.ratio[which(df.gm.lpa.ratio$gene_id %in% pc),"gene_type"] <- "mRNA"
  ggplot(df.gm.lpa.ratio, aes(y=RPKMSum.cyt,x=RPKMSum.nuc)) + geom_point()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_abline(slope=1,intercept=0)+
    ggtitle("K562 normailzed RPKM, cyt, nuc")
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratio-K562.pdf"))
  
  ggplot(df.gm.lpa.ratio, aes(y=RPKMSum.cyt,x=RPKMSum.nuc)) + geom_point()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_abline(slope=1,intercept=0) + 
    ggtitle("K562 normailzed RPKM, cyt, nuc") + xlim(0,10) + ylim(0,10)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratio-K562-zoom.pdf"))
  
  ggplot(df.gm.lpa.ratio, aes(log(RPKM.ratio))) + geom_density()  +
    facet_wrap(~gene_type,ncol=1) + theme_bw() + geom_vline(x=0) +
    ggtitle("K562 normailzed RPKM, cyt / nuc")  
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/cytNuc-ratioDensity-K562-.pdf"))
  
  df.pos <- df.together[which(df.together$transTotalRPKMnorm > 0 & df.together$intronTotalRPKMnorm),]
  df.pos$exp <- factor(paste0(df.pos$rnaExtract,".",df.pos$localization),
                       levels=c( "longNonPolyA.nucleus", "longPolyA.nucleus","longPolyA.cytosol","longNonPolyA.cytosol" ))
  
  ggplot(df.pos, aes(x=transTotalRPKMnorm,y=intronTotalRPKMnorm,color=factor(gene_type),fill=gene_type))+ geom_point()+
    theme_bw() + facet_wrap(rnaExtract ~ localization,scale="free") + 
    ggtitle("GM12878 + K562") +
    scale_color_manual(values=c("red","green","yellow"))
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/intron-v-trans-dotplot.pdf"))
  
  ggplot(df.pos[which(df.pos$cell == "GM12878"),],
         aes(x=((transTotalRPKMnorm )/(intronTotalRPKMnorm )),fill=gene_type))+ geom_density(alpha=I(0.4))+
    theme_bw() + facet_wrap(~exp,scale="free",ncol=1) + 
    ggtitle("GM12878  trans/intron normailzed RPKM") + xlim(0,40)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/intron-v-trans-density-GM12878.pdf"))
  
  ggplot(df.pos[which(df.pos$cell == "K562"),],
         aes(x=((transTotalRPKMnorm )/(intronTotalRPKMnorm )),fill=gene_type))+ geom_density(alpha=I(0.4))+
    theme_bw() + facet_wrap(~exp,scale="free",ncol=1) + 
    ggtitle("K562 trans/intron normailzed RPKM") + xlim(0,40)
  ggsave(getFullPath("plots/rnaExpr/mappedReads/starSpikeIn/intron-v-trans-density-K562.pdf"))
  
  
}


readSums <- function(df){
  df.mRNA = subset(df, gene_type == "mRNA")
  df.lncRNA = subset(df,gene_type == "lncRNA")
  df.other = subset(df, gene_type == "other")
  totalReads <- sum(df.mRNA$intronTotalReads,na.rm=TRUE) + sum(df.mRNA$transcriptTotalReads) +
    sum(df.lncRNA$transcriptTotalReads,na.rm=TRUE) + sum(df.lncRNA$intronTotalReads,na.rm=TRUE)
  
  
  l <- list(mRNA.ex.reads = sum(df.mRNA$transcriptTotalReads,na.rm=TRUE),
            mRNA.int.reads = sum(df.mRNA$intronTotalReads,na.rm=TRUE),
            mRNA.ex.readFrac = sum(df.mRNA$transcriptTotalReads)/totalReads,
            mRNA.int.readFrac = sum(df.mRNA$intronTotalReads,na.rm=TRUE)/totalReads,
            mRNA.ex.readPerLen = sum(df.mRNA$transcriptReadsPerLength,na.rm=TRUE),
            mRNA.int.readPerLen = sum(df.mRNA$intronReadsPerLength,na.rm=TRUE),
            mRNA.ex.RPKM =  sum(df.mRNA$transTotalRPKM,na.rm=TRUE),
            mRNA.int.RPKM = sum(df.mRNA$intronTotalRPKM,na.rm=TRUE),
            mRNA.ex.RPKMspikeIn =  sum(df.mRNA$transcriptTotalRPKM_spikeIn,na.rm=TRUE),
            mRNA.int.RPKMspikeIn = sum(df.mRNA$intronTotalRPKM_spikeIn,na.rm=TRUE),
            mRNA.ex.ReadsPerTrans =  sum(df.mRNA$transcriptTotalReads,na.rm=TRUE)/max(df.mRNA$transcriptCount,na.rm=TRUE),
            mRNA.int.ReadsPerTrans = sum(df.mRNA$intronTotalReads,na.rm=TRUE)/max(df.mRNA$intronCount,na.rm=TRUE),
            lncRNA.ex.reads = sum(df.lncRNA$transcriptTotalReads,na.rm=TRUE),
            lncRNA.int.reads = sum(df.lncRNA$intronTotalReads,na.rm=TRUE),
            lncRNA.ex.readFrac = sum(df.lncRNA$transcriptTotalReads,na.rm=TRUE)/totalReads,
            lncRNA.int.readFrac = sum(df.lncRNA$intronTotalReads,na.rm=TRUE)/totalReads,
            lncRNA.ex.readPerLen = sum(df.lncRNA$transcriptReadsPerLength,na.rm=TRUE),
            lncRNA.int.readPerLen = sum(df.lncRNA$intronReadsPerLength,na.rm=TRUE),     
            lncRNA.ex.RPKM =  sum(df.lncRNA$transTotalRPKM,na.rm=TRUE),
            lncRNA.int.RPKM = sum(df.lncRNA$intronTotalRPKM,na.rm=TRUE),
            lncRNA.ex.RPKMspikeIn =  sum(df.lncRNA$transcriptTotalRPKM_spikeIn,na.rm=TRUE),
            lncRNA.int.RPKMspikeIn = sum(df.lncRNA$intronTotalRPKM_spikeIn,na.rm=TRUE),
            lncRNA.ex.ReadsPerTrans =  sum(df.lncRNA$transcriptTotalReads,na.rm=TRUE)/max(df.lncRNA$transcriptCount,na.rm=TRUE),
            lncRNA.int.ReadsPerTrans = sum(df.lncRNA$intronTotalReads,na.rm=TRUE)/max(df.lncRNA$intronCount,na.rm=TRUE))
  
  nms <- names(l)
  nms.one <- t(sapply(nms, function(x)strsplit(x, "\\.")[[1]]))
  d <- cbind(as.data.frame(nms.one), as.numeric(l))
  colnames(d) <- c("transcript", "gene_regions", "measure","count")
  d
}



getReadCountsForCells <- function(){
  cells <- c("A549", "GM12878", "H1-hESC", "HeLa-S3", "HepG2", "HUVEC", "IMR90", "K562", "MCF-7", "NHEK", "SK-N-SH")
  sapply(cells, function(cell)readsCount(cell))
}


readsCount <- function(celltype="GM12878"){
  print(paste("celltype", celltype, "fetching data..."))
  df.together <- getTranscriptData(celltype=celltype, rnaExtract="longPolyA")
  df.together$replicate <- ifelse(df.together$replicate > 2, df.together$replicate - 2 , df.together$replicate)
  df.gm.cyt <- subset(df.together,cell==celltype &  localization=="cytosol" & rnaExtract == "longPolyA" )
  df.gm.nuc <- subset(df.together,cell==celltype &  localization =="nucleus" & rnaExtract == "longPolyA" )
  sum.gm.nuc.rep1 <- readSums(df.gm.nuc[which(df.gm.nuc$replicate == 1),])
  sum.gm.nuc.rep2 <- readSums(df.gm.nuc[which(df.gm.nuc$replicate == 2),])
  if((dim(df.gm.nuc[which(df.gm.nuc$replicate == 1),])[1] == 0 && (dim(df.gm.nuc[which(df.gm.nuc$replicate == 2),])[1] == 0))){
    stop("cannot get read counts for nuclear rnas-seq")
  }
  
  if(dim(df.gm.nuc[which(df.gm.nuc$replicate == 1),])[1] == 0 ){
    tmp <- df.gm.nuc
    tmp$replicate <- 1
    df.gm.nuc <- rbind(df.gm.nuc,tmp)
    rm(tmp)
    sum.gm.nuc.rep1 <- sum.gm.nuc.rep2
    sum.gm.nuc.rep1$tag = "N.A+.R2_copy"
    
    sum.gm.nuc.rep2$tag = "N.A+.R2"
  } else if(dim(df.gm.nuc[which(df.gm.nuc$replicate == 2),])[1] == 0 ){
    tmp <- df.gm.nuc
    tmp$replicate <- 2
    df.gm.nuc <- rbind(df.gm.nuc,tmp)
    rm(tmp)
    
    sum.gm.nuc.rep2 <- sum.gm.nuc.rep1
    sum.gm.nuc.rep2$tag = "N.A+.R1_copy"
    sum.gm.nuc.rep1$tag = "N.A+.R1"
    } else {
    sum.gm.nuc.rep1$tag = "N.A+.R1"
    sum.gm.nuc.rep2$tag = "N.A+.R2"
  }
  

 
  
  sum.gm.nuc.rep1$localization = "nucleus"
  sum.gm.nuc.rep2$localization = "nucleus"
  sum.gm.nuc.rep1$rnaExtract = "longPolyA"
  sum.gm.nuc.rep2$rnaExtract = "longPolyA"
  sum.gm.nuc.rep1$replicate = 1
  sum.gm.nuc.rep2$replicate =  2
  sum.gm.cyt.rep1 <- readSums(df.gm.cyt[which(df.gm.cyt$replicate == 1),])
  sum.gm.cyt.rep2 <- readSums(df.gm.cyt[which(df.gm.cyt$replicate == 2),])

  if((dim(df.gm.cyt[which(df.gm.cyt$replicate == 1),])[1] == 0 && (dim(df.gm.cyt[which(df.gm.cyt$replicate == 2),])[1] == 0))){
    stop("cannot get read counts for cytlear rnas-seq")
  }
  
  if(dim(df.gm.cyt[which(df.gm.cyt$replicate == 1),])[1] == 0 ){
    tmp <- df.gm.cyt
    tmp$replicate <- 1
    df.gm.cyt <- rbind(df.gm.cyt,tmp)
    rm(tmp)
    sum.gm.cyt.rep1 <- sum.gm.cyt.rep2
    sum.gm.cyt.rep1$tag = "C.A+.R2_copy"
    sum.gm.cyt.rep2$tag = "C.A+.R2"
  } else if(dim(df.gm.cyt[which(df.gm.cyt$replicate == 2),])[1] == 0 ){
    tmp <- df.gm.cyt
    tmp$replicate <- 2
    df.gm.cyt <- rbind(df.gm.cyt,tmp)
    rm(tmp)
    
    sum.gm.cyt.rep2 <- sum.gm.cyt.rep1
    sum.gm.cyt.rep2$tag = "C.A+.R1_copy"
    sum.gm.cyt.rep1$tag = "C.A+.R1"
  } else {
    sum.gm.cyt.rep1$tag = "C.A+.R1"
    sum.gm.cyt.rep2$tag = "C.A+.R2"
  }
  

  sum.gm.cyt.rep1$localization = "cytosol"
  sum.gm.cyt.rep2$localization = "cytosol"
  sum.gm.cyt.rep1$rnaExtract = "longPolyA"
  sum.gm.cyt.rep2$rnaExtract = "longPolyA"
  sum.gm.cyt.rep1$replicate = 1
  sum.gm.cyt.rep2$replicate =  2
  
  
  df.comb <- as.data.frame(rbind(sum.gm.nuc.rep1,sum.gm.nuc.rep2,sum.gm.cyt.rep1,sum.gm.cyt.rep2))
  df.comb$label = paste(df.comb$transcript,df.comb$gene_regions)
  
  # df.comb.reads <- df.comb[which(df.comb$measure == "reads"),]
  # df.comb.reads <- ddply(as.data.frame(df.comb.reads), .(tag), sumCount = sum(count))
  
  df.comb$labelPretty = c("mRNA.exon", "mRNA.intron", "lncRNA.exon", "lncRNA.intron")[as.numeric(sapply(df.comb$label,function(x) which(c("mRNA ex", "mRNA int", "lncRNA ex", "lncRNA int") %in% x)))]
  df.comb$region <- factor(df.comb$labelPretty, levels = rev(c("mRNA.exon", "mRNA.intron", "lncRNA.exon", "lncRNA.intron")))
  df.comb$symbolShort <- paste0(ifelse(df.comb$localization == "cytosol", "C","N"),".",ifelse(df.comb$rnaExtract == "longNonPolyA","A-","A+"))
  df.comb$symbol     <- paste0(  df.comb$symbolShort,".R",df.comb$replicate)
  symbolOrder <- c("N.A-","N.A+", "C.A+", "C.A-")
  df.comb$order1 <- sapply(df.comb$symbolShort, function(x)which(symbolOrder == x))
  df.comb$order2 <- df.comb$replicate
  
  
  plot.colors <- scale_fill_manual(values=c("mRNA.exon"= rgb(0.3765,0.5843,0.7882),
                                            "mRNA.intron"= rgb(0.4824,0.098,0.4745), 
                                            "lncRNA.exon"=rgb(0.3290,0.7094,0.7957),
                                            "lncRNA.intron"=rgb(0.9490,0.2667,0.4980)))
  
  
 ggplot(df.comb, aes(x=tag,y=count,fill=region)) + geom_bar(stat="bin") +
    facet_wrap(~measure,scale="free")+ plot.colors + theme_bw() + 
    scale_x_discrete(limits=unique(as.data.frame(arrange(df.comb,order1))[c("symbol")])$symbol) +
    ggtitle(paste0("celltype -> ", celltype,"\nFlux Capacitor mapped reads\ntotal reads = sum reads in all transcripts"))+
    thisTheme
  ggsave(paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"readCountsCytNuc-bars-",celltype,"-bars.pdf"),height=7,width=12)
  
  
  df.c <- as.data.frame(df.comb)
  df.c <- df.c[which(df.c$gene_regions == "ex"),]
  getExpectedRatio <- function(measure,transcript){
    nuc <- sum(df.c[which(df.c$measure == measure & df.c$localization == "nucleus" & df.c$transcript == transcript),"count"])
    cyt <- sum(df.c[which(df.c$measure == measure & df.c$localization == "cytosol" & df.c$transcript == transcript),"count"])
    log(cyt/nuc, base=exp(1))
  }
  getExpectedFrac <- function(measure,transcript){
    nuc <- sum(df.c[which(df.c$measure == measure & df.c$localization == "nucleus" & df.c$transcript == transcript),"count"])
    cyt <- sum(df.c[which(df.c$measure == measure & df.c$localization == "cytosol" & df.c$transcript == transcript),"count"])
    cyt/(cyt +nuc)
  }
  
  cyt.df <- merge(df.gm.cyt[which(df.gm.cyt$replicate == 1),],
                  df.gm.cyt[which(df.gm.cyt$replicate == 2),],
                  by="gene_id",suffixes= c(".rep1",".rep2"),all=TRUE)
  nuc.df <- merge(df.gm.nuc[which(df.gm.nuc$replicate == 1),],
                  df.gm.nuc[which(df.gm.nuc$replicate == 2),],
                  by="gene_id",suffixes= c(".rep1",".rep2"),all=TRUE)
  cCytNuc.df <- merge(cyt.df,nuc.df, by="gene_id", suffixes=c(".cyt",".nuc"))
  cCytNuc.df$nucExpressed <- 0
  cCytNuc.df$nucExpressed <- ifelse((cCytNuc.df$transcriptTotalReads.rep1.nuc + cCytNuc.df$transcriptTotalReads.rep2.nuc > 0), 1, 0)
  cCytNuc.df$cytExpressed <- ifelse((cCytNuc.df$transcriptTotalReads.rep1.cyt + cCytNuc.df$transcriptTotalReads.rep2.cyt > 0), 1, 0)
  
  
  
  
  cmRNA <- cCytNuc.df[which(cCytNuc.df$gene_type.rep2.nuc == "mRNA"),]
  venn.diagram(x=list(area1=which(cmRNA$cytExpressed > 0),
                      area2=which(cmRNA$nucExpressed > 0)),
               # cross.area=which( cCytNuc.df$nucExpressed +  cCytNuc.df$cytExpressed > 1 )),
               category=c("cytosol\nexpressed","nucleus\nexpressed"), 
               main.pos = c(0.5, 1.05),main=paste("celltype : ",celltype," mRNA"),
               height = 3100, width = 3000, resolution = 500,
               main.cex=2,
               filename=paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"vennDiagram-",celltype,"mRNA.tiff"))
  clncRNA <- cCytNuc.df[which(cCytNuc.df$gene_type.rep2.nuc == "lncRNA"),]
  venn.diagram(x=list(area1=which(clncRNA$cytExpressed > 0),
                      area2=which(clncRNA$nucExpressed > 0)),
               # cross.area=which( cCytNuc.df$nucExpressed +  cCytNuc.df$cytExpressed > 1 )),
               category=c("cytosol\nexpressed","nucleus\nexpressed"), 
               main.pos = c(0.5, 1.05),main=paste("celltype : ",celltype," lncRNA"),
               height = 3100, width = 3000, resolution = 500,
               main.cex=2,sub.cex=2,sub.col = "red",
               filename=paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"vennDiagram-",celltype,"lncRNA.tiff"))
  
  
  
  #dev.off()
  #graphics.off()
  
  ##cCytNucTransTotal.df <- cCytNuc.df[which(!(cCytNuc.df$transcriptTotalReads.rep1.cyt + cCytNuc.df$transcriptTotalReads.rep2.cyt > 0)),]
  #thrownOut.nuc <- sum(cCytNucTransTotal.df$transcriptTotalReads.rep1.nuc) + sum(cCytNucTransTotal.df$transcriptTotalReads.rep2.nuc)
  #cCytNucTransTotal.df <- cCytNuc.df[which(!(cCytNuc.df$transcriptTotalReads.rep1.nuc + cCytNuc.df$transcriptTotalReads.rep2.nuc > 0)),]
  #thrownOut.cyt <- sum(cCytNucTransTotal.df$transcriptTotalReads.rep1.cyt) + sum(cCytNucTransTotal.df$transcriptTotalReads.rep2.cyt)
  
  ##cytRepAdd.df <- ddply(cyt.df, .(gene_id),numcolwise(sum))
  #trans <- unique(cyt)
  
  cytNuc.df <- data.frame(gene_id = cyt.df$gene_id,
                          gene_type = cyt.df$gene_type.rep1,
                          readsCytFraction =  (cyt.df$transcriptTotalReads.rep1 + cyt.df$transcriptTotalReads.rep2)/(cyt.df$transcriptTotalReads.rep1 + cyt.df$transcriptTotalReads.rep2 + nuc.df$transcriptTotalReads.rep1 + nuc.df$transcriptTotalReads.rep2),
                          readsPerKbRatioCytFraction = (cyt.df$transcriptTotalReadsPerKb.rep1 + cyt.df$transcriptTotalReadsPerKb.rep2)/(cyt.df$transcriptTotalReadsPerKb.rep1 + cyt.df$transcriptTotalReadsPerKb.rep2 + nuc.df$transcriptTotalReadsPerKb.rep1 + nuc.df$transcriptTotalReadsPerKb.rep2),
                          RPKM_cytFraction =  (cyt.df$transTotalRPKM.rep1 + cyt.df$transTotalRPKM.rep2)/(cyt.df$transTotalRPKM.rep1 + cyt.df$transTotalRPKM.rep2 + nuc.df$transTotalRPKM.rep1 + nuc.df$transTotalRPKM.rep2),
                          RPKM_spikeIn_cytFraction = (cyt.df$transcriptTotalRPKM_spikeIn.rep1 + cyt.df$transcriptTotalRPKM_spikeIn.rep2)/ (cyt.df$transcriptTotalRPKM_spikeIn.rep1 + cyt.df$transcriptTotalRPKM_spikeIn.rep2 +nuc.df$transcriptTotalRPKM_spikeIn.rep1 + nuc.df$transcriptTotalRPKM_spikeIn.rep2))

  
  m.df <- melt(cytNuc.df, id=c("gene_id","gene_type"))
  m.df$expectedMean.mRNA <- 0
  m.df$expectedMean.lncRNA <- 0
  m.df[which(m.df$variable == "readsCytFraction"),"expectedMean.mRNA"] <- getExpectedFrac("reads","mRNA")
  m.df[which(m.df$variable == "readsCytFraction"),"expectedMean.lncRNA"] <- getExpectedFrac("reads","lncRNA")
  m.df[which(m.df$variable == "readsPerKbRatioCytFraction"),"expectedMean.mRNA"] <- getExpectedFrac("readPerLen","mRNA")
  m.df[which(m.df$variable == "readsPerKbRatioCytFraction"),"expectedMean.lncRNA"] <- getExpectedFrac("readPerLen","lncRNA")
  m.df[which(m.df$variable == "RPKM_cytFraction"),"expectedMean.mRNA"] <- getExpectedFrac("RPKM","mRNA")
  m.df[which(m.df$variable == "RPKM_cytFraction"),"expectedMean.lncRNA"] <- getExpectedFrac("RPKM","lncRNA")
  m.df[which(m.df$variable == "RPKM_spikeIn_cytFraction"),"expectedMean.mRNA"] <- getExpectedFrac("RPKMspikeIn","mRNA")
  m.df[which(m.df$variable == "RPKM_spikeIn_cytFraction"),"expectedMean.lncRNA"] <- getExpectedFrac("RPKMspikeIn","lncRNA")
  
  m.df$value <- ifelse(is.infinite(m.df$value),0,m.df$value)
  m.df$value <- ifelse(is.na(m.df$value),0,m.df$value)
  pos.df <-  subset(m.df,value > 0 & gene_type %in% c("lncRNA","mRNA"))
  ggplot(pos.df, aes(x=value,fill=gene_type))+ 
    geom_density(alpha=I(0.4))+ theme_bw() + facet_wrap(~variable)+
    thisTheme  + xlab("Frac = Cytosol / (Cytosol + Nucleus)") + xlim(-0.2,1.1) +
    ggtitle(paste("Celltype -> ",celltype,"\nFlux Capacitor mapped reads, combined replicates\nFrac = Cytosol / (Cytosol + Nucleus)"))
  ggsave(paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"fracCytNuc-density-",celltype,"-mappedReads.pdf"),height=7,width=6)
  
  #   geom_vline(data=pos.df, aes(xintercept=expectedMean.lncRNA,colour="green", linetype = "longdash"))+
  #    geom_vline(data=pos.df, aes(xintercept=expectedMean.mRNA,color="green")) 
  ggplot(pos.df, aes(x=value,fill=gene_type))+ 
    geom_bar(binwidth=0.05,position="dodge")+ theme_bw() + facet_wrap(~variable,scale="free")+
    thisTheme + xlab("Frac = Cytosol / (Cytosol + Nucleus)") + xlim(-0.2,1.1) +
    ggtitle(paste("Celltype -> ",celltype,"\nFlux Capacitor mapped reads, combined replicates\nFrac = Cytosol / (Cytosol + Nucleus)"))
  ggsave(paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"fracCytNuc-bars-",celltype,"-mappedReads.pdf"),height=7,width=6)
  
  readsFractionPseudo <- ecdf(cyt.df$transcriptTotalReads.rep1 + cyt.df$transcriptTotalReads.rep2 + nuc.df$transcriptTotalReads.rep1 + nuc.df$transcriptTotalReads.rep2)(0.05)
  readsPerKbRatioCytFractionPseudo  = ecdf(cyt.df$transcriptTotalReadsPerKb.rep1 + cyt.df$transcriptTotalReadsPerKb.rep2 + nuc.df$transcriptTotalReadsPerKb.rep1 + nuc.df$transcriptTotalReadsPerKb.rep2)(0.05)
  RPKM_cytFractionPseudo =  ecdf(cyt.df$transTotalRPKM.rep1 + cyt.df$transTotalRPKM.rep2 + nuc.df$transTotalRPKM.rep1 + nuc.df$transTotalRPKM.rep2)(0.05)
  RPKM_spikeIn_cytFractionPseudo = ecdf(cyt.df$transcriptTotalRPKM_spikeIn.rep1 + cyt.df$transcriptTotalRPKM_spikeIn.rep2 +nuc.df$transcriptTotalRPKM_spikeIn.rep1 + nuc.df$transcriptTotalRPKM_spikeIn.rep2)(0.05)
  cytNucPseudo.df <- data.frame(gene_id = cyt.df$gene_id,
                          gene_type = cyt.df$gene_type.rep1,
                          readsCytFraction =  (cyt.df$transcriptTotalReads.rep1 + cyt.df$transcriptTotalReads.rep2 + readsFractionPseudo)/(readsFractionPseudo +cyt.df$transcriptTotalReads.rep1 + cyt.df$transcriptTotalReads.rep2 + nuc.df$transcriptTotalReads.rep1 + nuc.df$transcriptTotalReads.rep2),
                          readsPerKbRatioCytFraction = (cyt.df$transcriptTotalReadsPerKb.rep1 + cyt.df$transcriptTotalReadsPerKb.rep2 + readsPerKbRatioCytFractionPseudo)/(readsPerKbRatioCytFractionPseudo + cyt.df$transcriptTotalReadsPerKb.rep1 + cyt.df$transcriptTotalReadsPerKb.rep2 + nuc.df$transcriptTotalReadsPerKb.rep1 + nuc.df$transcriptTotalReadsPerKb.rep2),
                          RPKM_cytFraction =  (cyt.df$transTotalRPKM.rep1 + cyt.df$transTotalRPKM.rep2 + RPKM_cytFractionPseudo)/(RPKM_cytFractionPseudo + cyt.df$transTotalRPKM.rep1 + cyt.df$transTotalRPKM.rep2 + nuc.df$transTotalRPKM.rep1 + nuc.df$transTotalRPKM.rep2),
                          RPKM_spikeIn_cytFraction = (cyt.df$transcriptTotalRPKM_spikeIn.rep1 + cyt.df$transcriptTotalRPKM_spikeIn.rep2 + RPKM_spikeIn_cytFractionPseudo)/ (RPKM_spikeIn_cytFractionPseudo + cyt.df$transcriptTotalRPKM_spikeIn.rep1 + cyt.df$transcriptTotalRPKM_spikeIn.rep2 +nuc.df$transcriptTotalRPKM_spikeIn.rep1 + nuc.df$transcriptTotalRPKM_spikeIn.rep2))
  mPseudo.df <- melt(cytNuc.df, id=c("gene_id","gene_type"))
  mPseudo.df$value <- ifelse(is.infinite(m.df$value),0,m.df$value)
  mPseudo.df$value <- ifelse(is.na(m.df$value),0,m.df$value)
  posPseudo.df <-  subset(m.df,value > 0 & gene_type %in% c("lncRNA","mRNA"))
  ggplot(posPseudo.df, aes(x=value,fill=gene_type))+ 
    geom_density(alpha=I(0.4))+ theme_bw() + facet_wrap(~variable)+
    thisTheme  + xlab("Frac = Cytosol / (Cytosol + Nucleus)") + xlim(-0.2,1.1) +
    ggtitle(paste("Celltype -> ",celltype,"\nFlux Capacitor mapped reads, combined replicates\nFrac = Cytosol + pseudo/ (pseudo +(Cytosol + Nucleus)\nwhere pseudo is 0.05 cdf of non-fraction value"))
  ggsave(paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"fracCytNucPseudo-density-",celltype,"-mappedReads.pdf"),height=7,width=6)
  
  #   geom_vline(data=pos.df, aes(xintercept=expectedMean.lncRNA,colour="green", linetype = "longdash"))+
  #    geom_vline(data=pos.df, aes(xintercept=expectedMean.mRNA,color="green")) 
  ggplot(posPseudo.df, aes(x=value,fill=gene_type))+ 
    geom_bar(binwidth=0.05,position="dodge")+ theme_bw() + facet_wrap(~variable,scale="free")+
    thisTheme + xlab("Frac = Cytosol / (Cytosol + Nucleus)") + xlim(-0.2,1.1) +
    ggtitle(paste("Celltype -> ",celltype,"\nFlux Capacitor mapped reads, combined replicates\nFrac = Cytosol + pseudo/ (pseudo +Cytosol + Nucleus)\nwhere pseudo is 0.05 cdf of non-fraction value"))
  ggsave(paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"fracCytNucPseudo-bars-",celltype,"-mappedReads.pdf"),height=7,width=6)
  
  summary.df <- group_by(pos.df, variable,gene_type) %.% 
                                 summarise(medianVal = median(value),
                                           meanVal   = mean(value),
                                           expectedmRNA =  mean(expectedMean.mRNA),
                                           expectedlncRNA = mean(expectedMean.lncRNA))
  
  summary.df <- as.data.frame(summary.df)
  summary.df$fracCytNucTotals <- ifelse(summary.df$gene_type == "lncRNA",summary.df$expectedlncRNA,
                                   summary.df$expectedmRNA)
  
  summary.df$expectedmRNA <- NULL;summary.df$expectedlncRNA <- NULL
  exportAsTable(file=paste0(getFullPath("plots/rnaExpr/mappedReads/starGrace/"),"fracCytNuc-summary-",celltype,".tab"),
                df=summary.df)
  
  #neg.df <-  subset(cyt.df,!(m.df$value > 0) & gene_type %in% c("lncRNA","mRNA"))
 # ddply(neg.df, .(variable,gene_type),summarise,sum(value))
  
}











