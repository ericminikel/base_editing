#!/bin/bash
#$ -l h_vmem=8G
#$ -cwd
#$ -q broad
#$ -l h_rt=48:00:00
#$ -j y
#$ -o qsub_out/redi.o$JOB_ID
#$ -M an@broadinstitute.org

# Load necessary modules
source /broad/software/scripts/useuse
reuse UGER
reuse Bowtie2
reuse Samtools

bowtie2-build /directory/amplicon.fasta /index_output_directory/index

bowtie2 -x /index_output_directory/index -U /directory/sample_name.fastq.gz -S /directory/sample_name.sam --very-fast-local

samtools view -S -b /directory/sample_name.sam > /directory/sample_name.bam

samtools sort /directory/sample_name.bam -o /directory/sample_name.sorted.bam

samtools index /directory/sample_name.sorted.bam

samtools idxstats /directory/sample_name.sorted.bam > /directory/sample_name.txt

input_bam="/directory/sample_name.sorted.bam"
reference_file="/directory/sample_name.txt"
output_folder="/output_directory"

# Create the output folder if it doesn't exist
mkdir -p ${output_folder}

reference_sequences=$(cat ${reference_file})

for ref in $reference_sequences; do
    output_bam="${output_folder}/${ref}.bam.gz"
    samtools view -h -b ${input_bam} ${ref} | gzip > ${output_bam}
done

input_folder="/directory/output_bam_files_sample_name"
output_folder="/directory/rhampseq-fastqs-sample_name"

# Create the output folder if it doesn't exist
mkdir -p ${output_folder}

# Process each .bam.gz file in the input folder
for bamfile in ${input_folder}/*.bam.gz; do
    # Extract the base name of the file (without path and extension)
    base_name=$(basename "${bamfile}" .bam.gz)

    # Define the output FASTQ file path
    fastq_file="${output_folder}/${base_name}.fastq"

    # Convert BAM to FASTQ
    echo "Processing ${bamfile} to ${fastq_file}"
    rm $fastq_file
    gunzip -c "${bamfile}" | samtools fastq >> "${fastq_file}"
done
