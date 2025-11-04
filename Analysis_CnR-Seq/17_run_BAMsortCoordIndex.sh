# script to sort BAM files by genomic coordinates and create index files

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_merged/*merged.bam`
do
	base=$(basename $sample ".bam")
	printf '%(%Y-%m-%d_%H:%M:%S)T\n' -1
	echo ${base} starting to sort...

	samtools sort ${dir}/mm10_merged/${base}.bam --threads 30 -o ${dir}/mm10_merged/${base}_sorted.bam

	echo ${base} finished sorting
done

cd ${dir}/mm10_merged/
samtools index -M *sorted.bam
