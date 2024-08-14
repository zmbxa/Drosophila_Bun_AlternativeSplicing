cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/data/clean

mkdir -p sam/summary
for i in Bun{2..4};do
echo $i
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 \
-x /storage/zhangyanxiaoLab/share/bowtie2_index/dm6 -1 fastq/${i}_R1.fq.gz -2 fastq/${i}_R2.fq.gz \
-S sam/${i}_bowtie2.sam &> sam/summary/${i}_bowtie2.txt
done

# or using ATAC pipeline?
conda activate paired_tag
/storage/zhangyanxiaoLab/share/Pipelines/atac-seq-pipeline-snakemake/bin/run_atacseq_vanilla.Zhanglab.sh -g dm6 -e niuyuxiao@westlake.edu -s 8


