rm(list=ls());gc()
setwd("~/projects/drosophila_CutTag_MaLab/analysis/")

# BiocManager::install("ChIPseeker")
library(ChIPseeker)
# BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
# BiocManager::install("org.Dm.eg.db")
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(ggplot2)
library(clusterProfiler)

Bun2_peak<-readPeakFile("macs2_peaks/Bun2_120_srt_narrow_peaks.narrowPeak") 
Bun2_peak = subset(Bun2_peak,nchar(as.character(seqnames))<6)

Bun3_peak<-readPeakFile("macs2_peaks/Bun3_120_srt_narrow_peaks.narrowPeak") 
Bun3_peak = subset(Bun3_peak,nchar(as.character(seqnames))<6)

Bun4_peak<-readPeakFile("macs2_peaks/Bun4_120_srt_narrow_peaks.narrowPeak") 
Bun4_peak = subset(Bun4_peak,nchar(as.character(seqnames))<6)

peaks = list(Bun2 = Bun2_peak,Bun3 = Bun3_peak,Bun4=Bun4_peak)
## visualize peaks
covplot(peaks)+ facet_grid(chr ~ .id)

# peak annotate
peakAnno <- annotatePeak(Bun2_peak,
                         tssRegion = c(-3000, 3000),  #启动子区域
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                         flankDistance = 5000,
                         #sameStrand = T,
                         TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                         annoDb = "org.Dm.eg.db") 
peakHeatmap(Bun2_peak,
            TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
            upstream=3000, downstream=3000,  #指定转录起始位点上下游
            color=rainbow(length(Bun2_peak)))
plotAvgProf2(Bun2_peak, 
             TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)

# peak visualization together
plotAvgProf2(peaks,TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
                          upstream=3000, downstream=3000,
                           xlab="Genomic Region (5'->3')", 
                           ylab = "Read Count Frequency",
                           conf = 0.95, resample = 1000)

## group annotate
peakAnnoList <- lapply(peaks, annotatePeak, 
                       TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                       tssRegion=c(-3000, 3000),
                       genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                       flankDistance = 5000,annoDb = "org.Dm.eg.db")
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

## annotated genes
genes <- lapply(peakAnnoList, function(i) 
  unique(as.data.frame(i)$geneId))
vennplot(genes,title="Venn of Nearest annotated genes")

# GO enrichment
# Run GO enrichment analysis 
entrez <- lapply(peakAnnoList, function(i) 
  unique(as.data.frame(i)$ENTREZID))
ego <- lapply(entrez, function(x){
  enrichGO(gene = x, 
           keyType = "ENTREZID", 
           OrgDb = org.Dm.eg.db, 
           ont = "BP", 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.05, 
           readable = TRUE)
})



## save results
dir.create("geneList",recursive = T)
for (i in names(peakAnnoList)) {
  print(i)
  write.table(as.data.frame(peakAnnoList[[i]]),col.names = T,row.names = F,quote = F,sep = "\t",file = paste0("geneList/",i,"_genelist.txt"))
  write.table(peakAnnoList[[i]]@annoStat,col.names = T,row.names = F,quote = F,sep = "\t",file = paste0("geneList/",i,"_annoStat.txt"))
}
write.table(Reduce(intersect,genes),file = "geneList/common_geneID.txt",col.names = F,row.names = F,sep = "\n",quote = F)

save.image("peak_stat.RData")
