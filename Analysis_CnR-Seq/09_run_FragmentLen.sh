# script to extract the 9th column of the SAM/BAM file which is the fragment length

for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_sorted/*_sorted.bam`
do
	dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"
	base=$(basename $sample "_sorted.bam")
	
	samtools view -F 0x04 ${sample} | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${dir}/fragmentLen/${base}_fragmentLen.txt

	echo ${base} finished
done
