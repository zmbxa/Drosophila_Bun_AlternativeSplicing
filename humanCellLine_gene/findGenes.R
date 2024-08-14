## read peak files and find binding genes
setwd("~/projects/drosophila_CutTag_MaLab/analysis/humanCellLine_gene/")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

peakList = list()
for(name in list.files("~/projects/drosophila_CutTag_MaLab/data/humanCellLine_gene/bed/")){
  peakList[gsub("\\.bed","",name)] = readPeakFile(paste0("~/projects/drosophila_CutTag_MaLab/data/humanCellLine_gene/bed/",name))
}

## define promoter length
promoter <- getPromoters(TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peakList, getTagMatrix, windows=promoter)
# peak annotation
peakAnnoList <- lapply(peakList, annotatePeak, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                       tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
# plot distribution
plotAnnoBar(peakAnnoList)

### save gene list
write.csv(as.data.frame(peakAnnoList[["TEAD1_HepG2_peakFromENCODE"]]),row.names = F,quote = T,file = paste0("gene_list/TEAD_HepG2_genelist.csv"))
write.csv(as.data.frame(peakAnnoList[["TEAD1_peakFromENCODE"]]),row.names = F,quote = T,file = paste0("gene_list/TEAD_genelist.csv"))
write.csv(as.data.frame(peakAnnoList[["TSC22D1_peakFromENCODE"]]),row.names = F,quote = T,file = paste0("gene_list/TSC22D1_genelist.csv"))






