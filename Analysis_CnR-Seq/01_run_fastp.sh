# script to run fastp trimming and filtering on all PE fastq files in a folder

cd ../raw_data
echo "changed directory to raw_data"

for first in `ls *_R1_001.fastq.gz`;
do
	second=`echo ${first} | sed 's/_R1/_R2/g'`
	prefix=`echo ${first} | sed 's/_R1_001.fastq.gz//g'`
	echo ${prefix} being processed ...
	fastp -i ${first} -I ${second} -o ${prefix}_R1_001_filtered.fastq.gz -O ${prefix}_R2_001_filtered.fastq.gz --trim_poly_g --low_complexity_filter --length_required 50 --qualified_quality_phred 30 --thread 20 --json ${prefix}_fastp.json --html ${prefix}_fastp.html
	echo ${prefix} finished
done

mkdir ../results/fastp_filtered_data
mkdir ../results/fastp_filtered_data/reports
mv *filtered.fastq.gz ../results/fastp_filtered_data
mv *.html ../results/fastp_filtered_data/reports
mv *.json ../results/fastp_filtered_data/reports
echo "results were moved to results/fastp_filtered_data"
