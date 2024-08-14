setwd("/storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/RNAseq")
rm(list = ls());gc()

counts = read.table("~/projects/drosophila_CutTag_MaLab/data/RNAseq/clean/te/te_counts.txt",header = T,row.names = 1)
colnames(counts) = limma::strsplit2(colnames(counts),"\\.")[,2]
counts=counts[,order(colnames(counts))]

y=DGEList(counts = counts,group = rep(2:1,each = 3))
y=y[which(rowSums(edgeR::cpm(y)>1)>=2),]
plotMDS(y,main="MDS plot of BUN and WT samples, RNAseq")

load("~/pipelines/RNAseq/analysis/get_edgeR_DE.RData")
BUNoverWT_out = get_EdgeR_DEout(counts,group = rep(2:1,each = 3),fdr = 0.05,abs_FC = 0)

ggplot(BUNoverWT_out,aes(x=logFC,y=-log10(FDR_tag),colour=change,label=gene_name))+
  geom_point(alpha=0.5, size=3.5)+scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggrepel::geom_text_repel(data=head(BUNoverWT_out,15),aes(label=gene_name),
                           size = 3,box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.8, "lines"), 
                           segment.color = "black", 
                           show.legend = FALSE)+
  ggtitle("BUN(3) over WT(3), RNA-seq")
enGO = compareCluster(gene_name~change,data = filter(BUNoverWT_out,change!="Stable"),
               fun = "enrichGO",keyType = "SYMBOL",OrgDb='org.Dm.eg.db',pvalueCutoff=0.05,ont = "BP")
enGO_Up = enrichGO(filter(BUNoverWT_out,change=="Up")$gene_name,OrgDb = "org.Dm.eg.db",keyType = "SYMBOL",
                   ont = "BP",pAdjustMethod = "BH")
barplot(enGO_Up,showCategory = 10)+ggtitle("Enriched GO pathway of BUN/WT Up-regulated genes")

# enrichGO(filter(BUNoverWT_out,change=="Down")$gene_name,OrgDb = "org.Dm.eg.db",keyType = "SYMBOL",
#          ont = "BP",pAdjustMethod = "none")

## saving results
write.csv(filter(BUNoverWT_out,change!="Stable"),file = "BUNoverWT_DEGenes.csv",quote = F,row.names = F)
write.csv(BUNoverWT_out,file = "BUNoverWT_out.csv",quote = F,row.names = F)
write.csv(enGO@compareClusterResult,file = "BUNoverWT_GOenrichment.csv",quote = F,row.names = F)
save.image("BUNoverWT_DEanalysis_RNAseq.RData")


