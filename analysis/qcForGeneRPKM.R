







getReportsForGene <- function(){
  list(flux=readInTable(getFullPath("/data/fluxCapData-lpa-proc-REPORT.tab")),
       rpkmFromBam=readInTable(getFullPath("/data/rpkmFromBAMCapData-lpa-proc-REPORT.tab")),
       rsem=readInTable(getFullPath("/data/rsemCapData-lpa-proc-REPORT.tab")))
}

plotRepQC <- function(){
  rep = getReportsForGene()
  colNames <- c("cell", "localization", "replicate","genesFound", "experiment")
  flux=rep$flux[colNames]
  rsem=rep$rsem[colNames]
  rpkmFromBam=rep$rpkmFromBam[colNames]
  flux$method = "flux"
  rsem$method = "rsem"
  rpkmFromBam$method = "rpkmFromBam"
  comb <- rbind(flux,rsem,rpkmFromBam)
  comb$tag <- with(comb,paste(cell,experiment,sep="."))

  ggplot(comb, aes(x=tag, y=genesFound))+facet_grid(method~.) + geom_bar() +
    ggtitle("Number of genes with read information\nFacets are different methods")+
    xlab("cell.localization.replicate") +
    ylab("count of genes found") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  ggsave(file=getFullPath("plots/rnaExpr/mappedReads/compareMethods/foundGenes.pdf"),height=7,width=12)


}













