# script to run bowtie2 on all filtered PE fastq files
# then convert SAM output to BAM and delete SAM

for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/fastp_filtered_data/*_R1_001_filtered.fastq.gz`
do
	dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/fastp_filtered_data"
	base=$(basename $sample "_R1_001_filtered.fastq.gz")
	printf '%(%Y-%m-%d_%H:%M:%S)T\n' -1
	echo ${base} starting alignment...
	bowtie2 -p 30 -q --end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 --dovetail -x /group/db/bowtie2_indexes/mm10 -1 ${dir}/${base}_R1_001_filtered.fastq.gz -2 ${dir}/${base}_R2_001_filtered.fastq.gz -S ${dir}/${base}.sam

	samtools view -h -S -b -o ${dir}/${base}.bam ${dir}/${base}.sam

	rm ${dir}/${base}.sam
	echo ${base} finished
done
