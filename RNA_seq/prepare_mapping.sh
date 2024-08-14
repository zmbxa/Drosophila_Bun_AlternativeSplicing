cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/data/RNAseq/clean/fastq
md5sum *   ##  grep fq.gz ../md5.txt 

## make rmsk gtf file
zcat ~/annotations/dm6.rmsk.txt.gz|awk '{OFS=FS="\t";a[$11]++;print$6,"dm6_rmsk","exon",$7,$8,$9,$10,".",$11,$11"_dup"a[$11]-1,$13,$12}'  | \
awk '{OFS=FS="\t";gsub("_dup0","",$10);print$1,$2,$3,$4,$5,$6,$7,$8,"gene_id ""\""$9"\"; ""transcript_id ""\""$10"\"; ""family_id ""\""$11"\"; ""class_id ""\""$12"\";"}' >/storage/zhangyanxiaoLab/niuyuxiao/annotations/gtf/dm6_rmsk.gtf

## rename for pipeline
rename '-' '' * -v # -n
rename fq fastq * -v # -n
rename _R _ * -v # -n
cd ..

conda activate snakemake
snakemake --cores 4 -s ~/pipelines/RNAseq/mapping/bulk_PE_RNA_dm6.reminder.smk --use-conda --conda-frontend conda --dryrun




