### stat TF co-binding result
setwd("~/projects/drosophila_CutTag_MaLab/analysis/ENCODE_ref/")
meta = read.delim("metadata.tsv")

ovl_coBind = read.delim("Bun234common.ovl.stat",header = F)
union_coBind = read.delim("Bun234union.ovl.stat",header = F)

colnames(ovl_coBind) = colnames(union_coBind) = c("length","nLine","accession")
ovl_coBind <- merge(ovl_coBind,dplyr::select(meta,File.accession ,Experiment.target),
                    by.x = "accession",by.y = "File.accession")
union_coBind <- merge(union_coBind,dplyr::select(meta,File.accession ,Experiment.target),
                    by.x = "accession",by.y = "File.accession")

# niuyuxiao@heterochromatin:~/projects/drosophila_CutTag_MaLab/analysis/ENCODE_ref$ awk '{sum+=($3-$2)}END{print sum}' Bun234_peak_union.bed
# 2552573
# niuyuxiao@heterochromatin:~/projects/drosophila_CutTag_MaLab/analysis/ENCODE_ref$ awk '{sum+=($3-$2)}END{print sum}' Bun234_peak_overlap.bed
# 1223191

ovl_coBind$percent = ovl_coBind$length/1223191*100
union_coBind$percent = union_coBind$length/2552573*100

### top 10 co-binding TF
ovl_coBind %>% arrange(desc(percent)) %>% head(10)
union_coBind %>% arrange(desc(percent)) %>% head(10)

write.table(ovl_coBind,file = "ovl_coBind.txt",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(union_coBind,file = "union_coBind.txt",col.names = F,row.names = F,quote = F,sep = "\t")

save.image("ENCODE_TF_stat.RData")
