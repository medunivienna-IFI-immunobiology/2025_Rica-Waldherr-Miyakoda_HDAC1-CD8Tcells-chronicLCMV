# script to run picard AddOrReplaceReadGroups on all sorted bam files

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_sorted/*_sorted.bam`
do
	base=$(basename $sample "_sorted.bam")
	java -jar /group/opt/picard/build/libs/picard.jar AddOrReplaceReadGroups -I ${dir}/mm10_sorted/${base}_sorted.bam -O ${dir}/mm10_sorted/${base}_sorted_Rgadded.bam -LB lib -PL Illumina -PU unit -SM ${base}
	echo ${base} finished
done
