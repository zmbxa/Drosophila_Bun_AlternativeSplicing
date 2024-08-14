## build TSS reference 
cd /storage/zhangyanxiaoLab/niuyuxiao/annotations

# pip install GetTss
GetTss -d ucsc -g /storage/zhangyanxiaoLab/share/gtf/dm6.ncbiRefSeq.gtf -t dm6.ncbiRefSeq_TSS.bed

cd ~/projects/drosophila_CutTag_MaLab/analysis

### generate bigWig using RPGC normalization
for i in Bun{2..4};do
echo $i
bamCoverage --bam bam/${i}.nodup.120.sorted.bam -o bigWig/bin25/${i}.nodup.120.sorted.25.bw --effectiveGenomeSize 129941135 --normalizeUsing RPGC --binSize 25 -p 10 -e 100
bamCoverage --bam bam/${i}.nodup.120.sorted.bam -o bigWig/bin10k/${i}.nodup.120.sorted.25.bw --effectiveGenomeSize 129941135 --normalizeUsing RPGC --binSize 10000 --smoothLength 100000 25 -p 10 -e 100
done

### compute matrix for 3 replicates
computeMatrix reference-point -R /storage/zhangyanxiaoLab/share/gtf/dm6.ncbiRefSeq.gtf  \
-S /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/bigWig/bin25/* \
-b 3000 -a 3000 --referencePoint TSS \
-p 8 --skipZeros --missingDataAsZero -o deeptools/all_TSS.mat.gz
 

computeMatrix scale-regions -R /storage/zhangyanxiaoLab/share/gtf/dm6.ncbiRefSeq.gtf \
-S /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/bigWig/bin25/* \
-b 3000 -a 3000 --regionBodyLength 5000 \
-p 8 --skipZeros --missingDataAsZero -o deeptools/all_region.mat.gz

### TSS ###
# plot profile
plotProfile -m deeptools/all_TSS.mat.gz -out deeptools/profileTSS.png --plotTitle "Bun2,3,4, profile" --perGroup --legendLocation upper-right --plotWidth 20 --plotHeight 15
# plot heatmap
plotHeatmap -m deeptools/all_TSS.mat.gz -out deeptools/heatmapTSS.png --samplesLabel Bun2 Bun3 Bun4 --heatmapWidth 3 --heatmapHeight 10
### region
plotProfile -m deeptools/all_region.mat.gz -out deeptools/profile_region.png --plotTitle "Bun2,3,4, profile" --perGroup --legendLocation upper-right --plotWidth 20 --plotHeight 15


plotHeatmap -m deeptools/all_TSS.mat.gz -out deeptools/heatmap_zoomTSS.png --samplesLabel Bun2 Bun3 Bun4 \
--heatmapWidth 6 --heatmapHeight 18 --legendLocation upper-right
plotHeatmap -m deeptools/all_TSS.mat.gz -out deeptools/heatmap_zoomTSS.svg --samplesLabel Bun2 Bun3 Bun4 \
--heatmapWidth 6 --heatmapHeight 18 --legendLocation upper-right --plotFileFormat svg --verbose



