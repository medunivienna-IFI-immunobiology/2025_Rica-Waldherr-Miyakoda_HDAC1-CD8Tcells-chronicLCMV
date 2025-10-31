#script to run picard MarkDuplicates on all sorted bam files

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_sorted/*_sorted.bam`
do
	base=$(basename $sample "_sorted.bam")
	java -jar /group/opt/picard/build/libs/picard.jar MarkDuplicates -I ${dir}/mm10_sorted/${base}_sorted.bam -O ${dir}/mm10_sorted/${base}_sorted_mDup.bam -M ${dir}/mm10_sorted/picard_summary/${base}_picard_mDup.txt
	echo ${base} finished
done
