#!/bin/bash
#$ -l h_vmem=8G          # Request 8 GB of memory
#$ -cwd                  # Run the job in the current working directory
#$ -q broad              # Submit to the broad queue
#$ -l h_rt=48:00:00      # Set a 48-hour runtime limit
#$ -j y                  # Merge standard error with standard output
#$ -o qsub_out/job_output.$JOB_ID.out  # Output file for job logs

# Load necessary modules
source /broad/software/scripts/useuse
reuse UGER
reuse Anaconda3

# Activate the Conda environment
source activate myenv

# Define variables for script parameters
SAMPLE_ID="sample_ID"
AMPLICON_GRNA_FILE="/directory/amplicon.tsv"
REF_BASE="C"
ALT_BASE="T"

# Run the Python script with the specified arguments
python /03_edit_eff.py $SAMPLE_ID $AMPLICON_GRNA_FILE $REF_BASE $ALT_BASE
