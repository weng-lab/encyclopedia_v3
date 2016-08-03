#Jill E. Moore
#Weng Lab
#UMass Medical School
#Updated August 2016

#ENCODE Encyclopedia Version 3
#Pipeline for generating enhancer-like regions from DNase peaks using DNase and H3K27ac signals

module load bedtools/2.25.0
module load python/2.7.5

d=$1

if [[ $2 == "mm9" ]]
then
TSS=~/Lab/Reference/Mouse/GencodeM1/TSS.4K.bed
ChromInfo=~/Lab/Reference/Mouse/mm9.ChromLength.txt
elif [[ $2 == "mm10" ]]
then
TSS=~/Lab/Reference/Mouse/GencodeM8/TSS.4K.bed
ChromInfo=~/Lab/Reference/Mouse/ChromInfo.txt
else
TSS=~/Lab/Reference/Human/Gencode19/TSS.4K.bed
ChromInfo=~/Lab/Reference/Human/chromInfo.txt
fi

echo $TSS $ChromInfo

dset=$(grep $d Data-Files | awk -F "\t" '{print $5}')
dsig=$(grep $d Data-Files | awk -F "\t" '{print $7}')
dpeaks=$(grep $d Data-Files | awk -F "\t" '{print $6}')

hset=$(grep $d Data-Files | awk -F "\t" '{print $1}')
hsig=$(grep $d Data-Files | awk -F "\t" '{print $3}')
hpeaks=$(grep $d Data-Files | awk -F "\t" '{print $2}')

echo /project/umw_zhiping_weng/0_metadata/encode/data/$dset/$dpeaks.bed.gz
cp /project/umw_zhiping_weng/0_metadata/encode/data/$dset/$dpeaks.bed.gz dpeaks.bed.gz
gunzip dpeaks.bed.gz

sort -rgk8,8 -rgk7,7 dpeaks.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" "Peak_"NR "\t" $5 "\t" $6 "\t" $7 "\t" $8}' > T

j=250
awk -F "\t" '{printf "%s\t%.0f\t%.0f\t%s\n", $1,($3+$2)/2-'$j',($3+$2)/2+'$j',$4}' T | awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4 ; else print $0}' | sort -u > little
~/bin/bigWigAverageOverBed -bedOut=out2.bed /project/umw_zhiping_weng/0_metadata/encode/data/$dset/$dsig.bigWig little out2

j=1000
awk -F "\t" '{printf "%s\t%.0f\t%.0f\t%s\n", $1,($3+$2)/2-'$j',($3+$2)/2+'$j',$4}' T | awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4 ; else print $0}' | sort -u > little
~/bin/bigWigAverageOverBed -bedOut=out3.bed /project/umw_zhiping_weng/0_metadata/encode/data/$hset/$hsig.bigWig little out3 

sort -k4,4 T > sorted1
sort -k1,1 out2 > sorted2
sort -k1,1 out3 > sorted3

paste sorted1 sorted2 sorted3 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $8 "\t" $13 "\t" $20}' > z

sort -rgk8,8 z | awk '{print $0 "\t" NR}' | sort -rgk9,9 | awk '{print $0 "\t" NR}' | awk '{print $0 "\t" ($10+$11)/2}' | sort -gk12,12  | awk '{print $0 "\t" NR}'> Ranked-List

if [[ $hpeaks == *"wgEncodeBroadHistone"* ]]
then
cp /project/umw_zhiping_weng/0_metadata/anshul/projects/encode/rawdata/peaks_histone/mar2012/narrow/combrep_and_ppr/$hpeaks ./hpeaks.bed.gz
gunzip hpeaks.bed.gz
else
cp /project/umw_zhiping_weng/0_metadata/encode/data/$hset/$hpeaks.bed.gz hpeaks.bed.gz
gunzip hpeaks.bed.gz
fi

sort -k1,1 -k2,2n hpeaks.bed > sorted4

bedtools merge -i sorted4 -d 250 > merge.bed
bedtools intersect -wo -f 1 -a Ranked-List -b merge.bed > First
bedtools intersect -v -f 1 -a Ranked-List -b merge.bed > no

bedtools merge -i sorted4 -d 500 > merge.bed
bedtools intersect -wo -f 1 -a no -b merge.bed > Second
bedtools intersect -v -f 1 -a no -b merge.bed > no2

bedtools intersect -wo -a no2 -b merge.bed > Third
bedtools intersect -v -a no2 -b merge.bed > no3

awk '{print $14 "\t" $15 "\t" $16 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13}' First > tmp
awk '{print $14 "\t" $15 "\t" $16 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13}' Second >> tmp
awk '{if ($2 < $15 && $3 > $16) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13; else if ($2 < $15) print $14 "\t" $2-500 "\t" $16 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13; else print $14 "\t" $15 "\t" $3+500 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13}' Third >> tmp
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" 1 "\t" "." "\t" $2 "\t" $3 "\t" $13}' no3 >> tmp

sort -k1,1 -k2,2n tmp > sorted-tmp
bedtools merge -i sorted-tmp -c 4 -o collapse > tmp-merge

python ~/Projects/ENCODE/Scripts/prune.dnase.enhancer.py tmp-merge sorted-tmp | sort -k9,9n | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > $d.EP

bedtools intersect -v -a $d.EP -b $TSS > EP
bedtools intersect -u -a $d.EP -b $TSS > P

awk '{print $1 "\t" $2 "\t" $3 "\t" "Distal-Prediction-"NR "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "255,205,0"}' EP | head -n 20000 > $dset'_'$hset'_predictions.bed'
m=$(head -n 20000 EP | tail -1 | awk '{print $9}')
awk '{if ($9 < '$m') print $1 "\t" $2 "\t" $3 "\t" "Proximal-Prediction-"NR "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "255,0,0"}' P >> $dset'_'$hset'_predictions.bed'

sort -k1,1 -k2,2n $dset'_'$hset'_predictions.bed' > tmp
~/bin/bedToBigBed tmp $ChromInfo $dset'_'$hset'_predictions.bigBed' -type=bed9
rm hpeaks.bed dpeaks.bed *tmp* no* First Second Third sorted* little out* T z Ranked-List merge.bed
