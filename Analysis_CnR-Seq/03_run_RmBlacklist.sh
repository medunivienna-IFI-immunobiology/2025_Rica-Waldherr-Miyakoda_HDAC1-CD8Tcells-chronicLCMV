# script to filter blacklisted regions

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls ${dir}/mm10/*.bam`
do
	base=$(basename $sample ".bam")
	bedtools intersect -v -abam ${dir}/mm10/${base}.bam -b ${dir}/mm10-blacklist.v2.bed > ${dir}/mm10_rmBL/${base}_RmBlacklist.bam
	echo ${base} finished
done
