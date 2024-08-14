cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis

mkdir bam bigWig
ln -s /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/data/clean/bam/* bam/
ln -s /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/data/clean/bigWig/* bigWig/

# sorting bam
for i in Bun{2..4};do
echo $i
sambamba sort -p -t 6 bam/${i}.nodup.bam 
done

# prepare bedgraph for peak calling
mkdir bdg
for i in Bun{2..4};do
bedtools genomecov -bg -ibam bam/${i}.nodup.sorted.bam > bdg/${i}.nodup.sorted.bedgraph
done

## SEACR call peak
for i in Bun{2..4};do
SEACR_1.3.sh bdg/${i}.nodup.sorted.bedgraph 0.05 norm relaxed peaks/$i
SEACR_1.3.sh bdg/${i}.nodup.sorted.bedgraph 0.05 norm stringent peaks/$i
done 


### follow https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10344426/#s4title
# keep reads below 120bp...?
for i in Bun{2..4};do
echo $i
samtools view -h bam/${i}.nodup.sorted.bam | LC_ALL=C awk -f /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/script/filter120.awk |samtools view -Sb - > bam/${i}.nodup.120.bam
sambamba index bam/${i}.nodup.120.bam -t 6 
done

## shift for cuttag(? ### donot shift!!
conda activate deeptools
for i in Bun{2..4};do
echo $i
alignmentSieve --numberOfProcessors 8 --ATACshift --bam bam/${i}.nodup.120.bam -o bam/${i}.nodup.120.shift.tmp.bam
# sort bam again
sambamba sort -t 8 -p bam/${i}.nodup.120.shift.tmp.bam -o bam/${i}.nodup.120.shift.bam
rm bam/${i}.nodup.120.shift.tmp.bam
done

## MACS2 call peak
for i in Bun{2..4};do
echo $i
#macs2 callpeak -t bam/${i}.nodup.120.shift.bam -g 1.2e8 -f BAMPE -n ${i}_120_shift_narrow --outdir macs2_peaks -q 0.01 -B --SPMR --keep-dup all 2> macs2_peaks/${i}.120.shift.macs2.narrow.log
#macs2 callpeak -t bam/${i}.nodup.120.shift.bam -g 1.2e8 -f BAMPE -n ${i}_120_shift_broad --outdir macs2_peaks -q 0.01 --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> macs2_peaks/${i}.120.shift.macs2.broad.log
macs2 callpeak -t bam/${i}.nodup.sorted.bam -g 1.2e8 -f BAMPE -n ${i}_narrow --outdir macs2_peaks -q 0.01 -B --SPMR --keep-dup all 2> macs2_peaks/${i}.macs2.narrow.log
macs2 callpeak -t bam/${i}.nodup.sorted.bam -g 1.2e8 -f BAMPE -n ${i}_broad --outdir macs2_peaks -q 0.01 --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all 2> macs2_peaks/${i}.macs2.broad.log
done 

## donnot need shift because its paired-end data; reserve no-shift narrow results
# but 120bp filter?
for i in Bun{2..4};do
echo $i
sambamba sort bam/${i}.nodup.120.bam -p -t 8 
macs2 callpeak -t bam/${i}.nodup.120.sorted.bam -g 1.2e8 -f BAMPE -n ${i}_120_srt_narrow --outdir macs2_peaks -q 0.01 -B --SPMR --keep-dup all 2> macs2_peaks/${i}_120_srt.macs2.narrow.log
done







