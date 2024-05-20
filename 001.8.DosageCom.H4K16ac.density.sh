#G4 density
##1.make 1kb windows
cd /home/yuss/flyG4/result/PQS/
bedtools makewindows -g /home/yuss/flyG4/data/ref/dm6.genome -w 1000 | bedtools intersect -a - -b /home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed -wa -c > 001.8.window1kb.MergedG4.bedGraph
bedtools makewindows -g /home/yuss/flyG4/data/ref/dm6.genome -w 1000 | bedtools intersect -a - -b 001.6.dm6.K.bed -wa -c > 001.8.window1kb.K.bedGraph
bedtools makewindows -g /home/yuss/flyG4/data/ref/dm6.genome -w 1000 | bedtools intersect -a - -b 001.6.dm6.PDS.bed -wa -c > 001.8.window1kb.PDS.bedGraph
##2.转成bw文件
for i in MergedG4 K PDS; do /home/qians/Biosoft/bedGraphToBigWig 001.8.window1kb.$i.bedGraph /home/yuss/flyG4/data/ref/dm6.genome 001.8.window1kb.$i.bw;done
##3.X A peak file
for i in male female; do grep -E "2L|2R|3L|3R" /home/yuss/flyG4/result/PQS/001.7.${i}_W_peaks.bed > 001.8.${i}.A.peaks.bed; done
for i in male female; do grep -E "X" /home/yuss/flyG4/result/PQS/001.7.${i}_W_peaks.bed > 001.8.${i}.X.peaks.bed; done
##4.computmatrix   
for i in MergedG4 K PDS; do for j in male female; do computeMatrix scale-regions -S 001.8.window1kb.${i}.bw -R 001.8.${j}.X.peaks.bed 001.8.${j}.A.peaks.bed --beforeRegionStartLength 500 --regionBodyLength 1000 --afterRegionStartLength 500 --skipZeros -o 001.8.${i}.${j}.matrix.gz; done; done
for i in MergedG4 K PDS; do for j in male female; do plotProfile -m 001.8.${i}.${j}.matrix.gz -out /home/yuss/flyG4/result/PQS/Picture/001.8.${i}.${j}.density.pdf --startLabel "Start" --endLabel "End"  --colors "#d6604d" "#06577A" --plotHeight 10 --plotWidth 7.5  --yAxisLabel "G4 Density" --samplesLabel $j --plotFileFormat pdf; done; done
##5.exmaple
computeMatrix scale-regions -S 001.8.window1kb.MergedG4.bw -R 001.8.male.X.peaks.bed --beforeRegionStartLength 500 --regionBodyLength 1000 --afterRegionStartLength 500 --skipZeros -o 001.8.MergedG4.male.X.matrix.gz
plotProfile -m 001.8.MergedG4.male.X.matrix.gz --startLabel "Start" --endLabel "End" --colors "#d6604d" "#06577A" --plotHeight 10 --plotWidth 8 --yAxisLabel "Density" --samplesLabel "Male" -out /home/yuss/flyG4/result/PQS/Picture/001.8.MergedG4.male.X.matrix.pdf --plotFileFormat pdf
