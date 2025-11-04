# script to sort BAM files by genomic coordinates and create index files

dir="/home/mwaldherr/Paper_RicaWaldherr_2025/analysis_cutnrun/results/bowtie2"

for sample in `ls $dir/mm10_merged/*merged.sorted.bam`
do
	base=$(basename $sample ".bam")

    bamCoverage -b ${dir}/mm10_merged/${base}.bam -o ${dir}/mm10_merged/${base}.bedgraph -of bedgraph -bs 20 --normalizeUsing None -p 30

	bamCoverage -b ${dir}/mm10_merged/${base}.bam -o ${dir}/mm10_merged/${base}.rpkm.BigWig -of bigwig -bs 20 --normalizeUsing RPKM -p 30

done
