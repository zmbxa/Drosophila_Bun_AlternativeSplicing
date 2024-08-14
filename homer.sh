## configuration
cd ~/tools/homer
perl configureHomer.pl -list|grep fly
perl configureHomer.pl -install 
perl configureHomer.pl -install dm6
perl configureHomer.pl -install fly-p

## prepare peak file from masc2
cd /storage/zhangyanxiaoLab/niuyuxiao/projects/drosophila_CutTag_MaLab/analysis/homer
for i in Bun{2..4};do 
echo $i
awk -v OFS="\t" '{print $4,$1,$2,$3,$6}' ../macs2_peaks/${i}_120_srt_narrow_peaks.narrowPeak >  peak_macs2/${i}.peaks
done

## enrich motif or co-factor?
# try Bun4 first
findMotifsGenome.pl  peak_macs2/Bun4.peaks dm6 bun4/ -size 120 -mask

###  motif enrichment for each replicates 
for i in Bun{2..3};do 
echo $i
mkdir -p $i
findMotifsGenome.pl  peak_macs2/${i}.peaks dm6 ${i}/ -size 120 -mask
done
# and common/all peaks of 3 replicates
cd ~/projects/drosophila_CutTag_MaLab/analysis/
awk -v OFS="\t" '{print "Overlap_peak"NR,$1,$2,$3,$6}' ENCODE_ref/Bun234_peak_overlap.bed >  homer/peak_macs2/Bun234ovl.peaks
awk -v OFS="\t" '{print "Union_peak"NR,$1,$2,$3,"."}' ENCODE_ref/Bun234_peak_union.bed > homer/peak_macs2/Bun234union.peaks
for i in Bun234*peaks;do 
name=$(basename $i)
name=${type%%p*}
echo $name
mkdir -p $i
findMotifsGenome.pl  peak_macs2/${i}.peaks dm6 ${i}/ -size 120 -mask
done


