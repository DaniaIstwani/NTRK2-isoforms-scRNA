#!/bin/bash
#SBATCH --job-name=dropest_job_both 
#SBATCH --output=dropest_output_both.txt
#SBATCH --error=dropest_error_both.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=bigmem

#SBATCH --array=0-31


#load modules:

# Variables:
samples=("10X52_1.bam.1" "10X22_2.bam.1" "10X28_3.bam.1" "10X28_2.bam.1" "10X50_4.bam.1" "10X52_3.bam.1" "10X50_1.bam.1" "10X19_2.bam.1" "10X05_1.bam.1" "10X52_4.bam.1" "10X49_1.bam.1" "10X52_2.bam.1" "10X05_2.bam.1" "10X50_3.bam.1"  
"10X22_1.bam.1" "10X36_3.bam.1" "10X20_2.bam.1" "10X24_2.bam.1" "10X06_2.bam.1" "10X35_2.bam.1" "10X20_1.bam.1" "10X36_2.bam.1" "10X36_1.bam.1" "10X38_2.bam.1" "10X35_1.bam.1" "10X06_1.bam.1" "10X07_1.bam.1" "10X38_1.bam.1"  
"10X38_3.bam.1" "10X87_2.bam.1" "10X87_1.bam.1" "10X86_3.bam.1")

mf_config="~/.Arcitecta/mflux.cfg"

sample_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/samples"

results_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/count_matrices/cm_both"

output_file="/data/gpfs/projects/punim2183/samples_processing/count_matrix_files/both/count_matrix.rds"

gtf_file="/data/gpfs/projects/punim2183/samples_processing/both.gtf"


# loop over each sample
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"

    module purge

    module load Java/17.0.6
    module load unimelb-mf-clients
    
    #fetch file from Mdediaflux:
    unimelb-mf-download --mf.config ~/.Arcitecta/mflux.cfg --out $(pwd) "/projects/proj-6030_ntrk2isoforms-1128.4.1092/samples/$sample"
    module purge

    sleep 5

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

done

# Log job's resource usage stats
JOBID=${SLURM_JOB_ID}

if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
    JOBID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
fi

my-job-stats -a -n -j "$JOBID"

