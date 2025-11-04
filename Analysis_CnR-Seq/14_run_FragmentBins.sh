# script to infer to which 500bp bin a fragment belongs (use midpoint of fragment)

binLen=500 
dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final"

for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final/*fragments.bed`
do
	base=$(basename $sample ".fragments.bed")
	awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${dir}/${base}.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${dir}/${base}.fragmentsCount.bin$binLen.bed
	
	echo ${base} finished
done
