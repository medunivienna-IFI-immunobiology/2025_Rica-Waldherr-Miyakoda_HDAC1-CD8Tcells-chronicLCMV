# script to convert mapped and sorted bam into bed file format
# additionally get clean bed file with only read pairs on same chromosome and fragment length <1000bp
# plus bed file containing only fragment related columns

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final"

for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final/*_sorted.bam`
do
	base=$(basename $sample ".bam")

	# convert into bed file format
	bedtools bamtobed -i ${dir}/${base}.bam -bedpe >${dir}/${base}.bed

	# keep the read pairs that are on the same chromosome and fragment length is <1000bp
	awk '$1==$4 && $6-$2 < 1000 {print $0}' ${dir}/${base}.bed >${dir}/${base}.clean.bed

	# only extract the fragment related columns
	cut -f 1,2,6 ${dir}/${base}.clean.bed | sort -k1,1 -k2,2n -k3,3n >${dir}/${base}.fragments.bed

	echo ${base} finished
done
