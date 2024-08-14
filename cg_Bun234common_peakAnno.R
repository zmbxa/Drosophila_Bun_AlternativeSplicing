### find shared peak of 3 replicates common with TF cg-dmelanogaster (ENCFF619NUF) binding genes
rm(list = ls());gc()
setwd("~/projects/drosophila_CutTag_MaLab/analysis/")

library(ChIPseeker)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)


cg_ovl = readPeakFile("ENCODE_ref/ovl_coBind/ENCFF619NUF_Bun234common.ovl.bed")

# peak annotate
peakAnno <- annotatePeak(cg_ovl,
                         tssRegion = c(-3000, 3000),  #启动子区域
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"),
                         flankDistance = 5000,
                         #sameStrand = T,
                         TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                         annoDb = "org.Dm.eg.db") 
peakHeatmap(cg_ovl,
            TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
            upstream=3000, downstream=3000,  #指定转录起始位点上下游
            color=rainbow(length(cg_ovl)))
plotAvgProf2(cg_ovl, 
             TxDb=TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
             upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", 
             ylab = "Read Count Frequency",
             conf = 0.95, resample = 1000)

plotAnnoBar(peakAnno)+ggtitle("Intersect peak of cg and 3 replicates common")
dir.create("TF_cg_ovl")
ggsave(filename = "TF_cg_ovl/cg_ovl_percentBar.jpeg",width = 250,height = 100,units = "mm")
## annotated genes
as.data.frame(peakAnno)

entrez <- unique(as.data.frame(peakAnno)$ENTREZID)
ego <- enrichGO(gene = entrez, 
                keyType = "ENTREZID", 
                OrgDb = org.Dm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
barplot(ego, showCategory = 20,font.size = 15)+ggtitle("Binding gene GO enrichment - (cg and 3 common)",)
ggsave(filename = "TF_cg_ovl/cg_ovl_geneEnrich.jpeg",width = 250,height = 500,units = "mm")
write.table(as.data.frame(peakAnno),col.names = T,row.names = F,quote = F,sep = "\t",file = paste0("TF_cg_ovl/cg_Bun234common_genelist.txt"))

