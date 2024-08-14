cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/ENCODE_ref
aria2c -c -x 4 -s 4 --max-concurrent-downloads=5 -i DownloadURL.txt  -d peak_bed

### find 3 replicates union and overlap
# union 
cd ..
cat macs2_peaks/Bun2_120_srt_narrow_peaks.narrowPeak macs2_peaks/Bun3_120_srt_narrow_peaks.narrowPeak macs2_peaks/Bun4_120_srt_narrow_peaks.narrowPeak |sort -k1,1 -k2,2n -|bedtools merge -i - >ENCODE_ref/Bun234_peak_union.bed
# overlap
bedtools intersect -a macs2_peaks/Bun2_120_srt_narrow_peaks.narrowPeak -b macs2_peaks/Bun3_120_srt_narrow_peaks.narrowPeak |bedtools intersect -a - -b macs2_peaks/Bun4_120_srt_narrow_peaks.narrowPeak >ENCODE_ref/Bun234_peak_overlap.bed

cd ENCODE_ref
mkdir -p ovl_coBind union_coBind
for i in peak_bed/*;do
name=$(basename $i)
name=${name%%.b*}
echo $name
bedtools intersect -a $i -b Bun234_peak_overlap.bed -wo > ovl_coBind/${name}_Bun234common.ovl.bed
bedtools intersect -a $i -b Bun234_peak_union.bed -wo > union_coBind/${name}_Bun234union.ovl.bed
done

for i in peak_bed/*;do
name=$(basename $i)
name=${name%%.b*}
echo $name
awk '{sum+=$21}END{print sum, NR}' ovl_coBind/${name}_Bun234common.ovl.bed >>Bun234common.ovl.stat
awk '{sum+=$14}END{print sum, NR}' union_coBind/${name}_Bun234union.ovl.bed >>Bun234union.ovl.stat
done

sed -i "s/ /\t/" Bun234union.ovl.stat
ls peak_bed/ | paste Bun234union.ovl.stat - > union.ovl.stat
sed -i "s/ /\t/" Bun234common.ovl.stat
ls peak_bed/ | paste Bun234common.ovl.stat - > common.ovl.stat

sed -i "s/.bed.gz//" Bun234union.ovl.stat
sed -i "s/.bed.gz//" Bun234common.ovl.stat

mv union.ovl.stat Bun234union.ovl.stat
mv common.ovl.stat Bun234common.ovl.stat
