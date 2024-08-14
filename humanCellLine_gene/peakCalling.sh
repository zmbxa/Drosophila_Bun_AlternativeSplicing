## convert bigWig to bed
###
cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/humanCellLine_gene
ln -s   /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/data/humanCellLine_gene/bam ./

### peak calling
mkdir macs2_out
for file in bam/*;do
name=${file%%.ba*}
name=$(basename $name)
echo $name
macs2 callpeak -t $file -f BAM -g hs  -n $name --outdir macs2_out --SPMR  --min-length 100  2> ${name}_macs.log
done

### use bed from ENCODE
pigz





### find common region between replicates
for type in TEAD1_HepG2 TEAD1 TSC22D1;do
echo $type
bedtools intersect -a macs2_out/${type}_rep1.narrowPeak -b macs2_out/${type}_rep2.narrowPeak 

