# script to merge all replicates per sample

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final"

newlist=""
for sample in `ls ~/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2/mm10_final/*KOctrl.mapped.bam`

do
	base=$(basename $sample "_RmBlacklistKOctrl.mapped.bam")
	base2=${base%-rep*}
	newlist="${newlist} ${base2}"
done
merged=`echo $newlist | tr " " "\n" | sort -u | tr "\n" " "`
echo $merged

for mergedsample in $merged
do
	## merge replicates to one merged.bam
	samtools merge -o ${dir}/${mergedsample}.KOctrl.merged.bam -n ${dir}/${mergedsample}-rep1_RmBlacklistKOctrl.mapped.bam ${dir}/${mergedsample}-rep2_RmBlacklistKOctrl.mapped.bam ${dir}/${mergedsample}-rep3_RmBlacklistKOctrl.mapped.bam

	echo ${mergedsample} finished
done
