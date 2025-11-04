# script to sort BAM files by names

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_sorted/*_sorted_mDup.bam`
do
	base=$(basename $sample "_sorted_mDup.bam")
	printf '%(%Y-%m-%d_%H:%M:%S)T\n' -1
	echo ${base} starting to sort...

	samtools sort -n --threads 20 ${dir}/mm10_sorted/${base}_sorted_mDup.bam  -o ${dir}/mm10_final/${base}.bam

	echo ${base} finished
done

# script to sort BAM files by names for controls

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_sorted/KO_controls/*_sorted_rmDup.bam`
do
        base=$(basename $sample "_sorted_rmDup.bam")
        printf '%(%Y-%m-%d_%H:%M:%S)T\n' -1
        echo ${base} starting to sort...

        samtools sort -n --threads 20 ${dir}/mm10_sorted/KO_controls/${base}_sorted_rmDup.bam  -o ${dir}/mm10_final/${base}KOctrl.bam

        echo ${base} finished
done

# script to keep only mapped and proper paired reads and header
for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final/*KOctrl.bam`
do
        base=$(basename $sample ".bam")

        ## filter and keep the mapped read pairs
        samtools view --threads 20 -h -bS -F 0x04 -f 0x02 ${sample} >${dir}/${base}.mapped.bam

        echo ${base} finished
done

