#!/bin/bash
#SBATCH --job-name=dropest_job_orig 
#SBATCH --output=dropest_output_orig.txt
#SBATCH --error=dropest_error_orig.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --partition=bigmem

#load modules:

# Variables

samples=("10X07_1.bam.1")

mf_config="~/.Arcitecta/mflux.cfg"

sample_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/sample"

results_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/count_matrices/cm_orig"

output_file="/data/gpfs/projects/punim2183/samples_processing/count_matrix_files/orig/count_matrix.rds"

gtf_file="/data/gpfs/projects/punim2183/samples_processing/orig.gtf"



echo "Processing sample: $sample"

module purge

module load Java/17.0.6
module load unimelb-mf-clients
    
    #fetch file from Mdediaflux:
unimelb-mf-download --mf.config ~/.Arcitecta/mflux.cfg --out $(pwd) "/projects/proj-6030_ntrk2isoforms-1128.4.1092/samples/$sample"
module purge


module load Java/11.0.18
module load foss/2022a
module load dropEst/0.8.6

    #run dropEst:
dropest -f -g "$gtf_file" -L e -c config.xml -o "$output_file" "$sample"
module purge



    # Rename the count matrix file
output_file_renamed="${output_file}_${sample}_$(basename $gtf_file).rds"
mv "$output_file" "$output_file_renamed"

   # Upload to MediaFlux
module load Java/17.0.6
module load unimelb-mf-clients

unimelb-mf-upload --mf.config "$mf_config" --csum-check --dest "$results_path" "$output_file_renamed"

  # Check if upload was successful before deleting the sample
if [ $? -eq 0 ]; then
    echo "Upload successful. Deleting sample: $sample"
    rm -rf "$(pwd)/$sample"
else
    echo "Upload failed. Sample not deleted: $sample"
fi

echo "Finished processing sample: $sample"

# Log job's resource usage stats
JOBID=${SLURM_JOB_ID}

if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
    JOBID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi

my-job-stats -a -n -j "$JOBID"

