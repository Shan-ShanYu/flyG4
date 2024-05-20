#Structure
cd /home/yuss/flyG4/result/PQS/
##1. get bedgraph
awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$5}' /home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed | bedtools sort -i - | bedtools merge -i - -c 4 -o max > 001.9.MergedG4.Structure.bedGraph
bedtools sort -i 001.6.dm6.K.bed | bedtools merge -i - -c 4 -o max > 001.9.K.Structure.bedGraph
bedtools sort -i 001.6.dm6.PDS.bed | bedtools merge -i - -c 4 -o max > 001.9.PDS.Structure.bedGraph
##2.to bw
/home/qians/Biosoft/bedGraphToBigWig 001.9.MergedG4.Structure.bedGraph /home/yuss/flyG4/data/ref/dm6.genome 001.9.MergedG4.Structure.bw
/home/qians/Biosoft/bedGraphToBigWig 001.9.K.Structure.bedGraph /home/yuss/flyG4/data/ref/dm6.genome 001.9.K.Structure.bw
/home/qians/Biosoft/bedGraphToBigWig 001.9.PDS.Structure.bedGraph /home/yuss/flyG4/data/ref/dm6.genome 001.9.PDS.Structure.bw
##3.computmatrix
for i in MergedG4 K PDS; do for j in male female; do computeMatrix scale-regions -S 001.9.${i}.Structure.bw -R 001.8.${j}.X.peaks.bed 001.8.${j}.A.peaks.bed --beforeRegionStartLength 500 --regionBodyLength 1000 --afterRegionStartLength 500 --skipZeros --binSize 20 -o 001.9.${i}.${j}.strcuture.matrix.gz; done; done
for i in MergedG4 K PDS; do for j in male female; do plotProfile -m 001.9.${i}.${j}.strcuture.matrix.gz -out /home/yuss/flyG4/result/PQS/Picture/001.9.${i}.${j}.strcuture.pdf --startLabel "Start" --endLabel "End"  --colors "#d6604d" "#06577A"  --plotHeight 10 --plotWidth 7.5  --yAxisLabel "Structural stability" --samplesLabel $j --plotFileFormat pdf; done; done
