#script to run picard MarkDuplicates on all sorted bam files and remove duplicates

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_sorted/KO_controls/*_sorted.bam`
do
	base=$(basename $sample "_sorted.bam")
	java -jar /group/opt/picard/build/libs/picard.jar MarkDuplicates -I ${dir}/mm10_sorted/KO_controls/${base}_sorted.bam -O ${dir}/mm10_sorted/KO_controls/${base}_sorted_rmDup.bam -M ${dir}/mm10_sorted/picard_summary/${base}_picard_rmDup.txt --REMOVE_DUPLICATES true
	echo ${base} finished
done
