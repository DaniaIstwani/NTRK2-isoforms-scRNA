#!/bin/bash
#SBATCH --job-name=dropest_job 
#SBATCH --output=dropest_output.txt
#SBATCH --error=dropest_error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --partition=bigmem


#load modules:

# Variables

sample="10X51_4.bam.1"

mf_config="~/.Arcitecta/mflux.cfg"

sample_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/samples"

results_path="/projects/proj-6030_ntrk2isoforms-1128.4.1092/count_matrices"

output_file="$(pwd)/count_matrix.rds"

gtf_file="$(pwd)/orig_s.gtf"

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
dropest -f -g "$gtf_file" -L e -c config.xml -o "$output_file" "$(pwd)/$sample"
module purge


#rename the count matrix rds
output_file_renamed="${output_file}_${sample}_$(basename $gtf_file).rds"
mv "$output_file" "$output_file_renamed"


module load Java/17.0.6
module load unimelb-mf-clients
#upload to Mf: 
unimelb-mf-upload --mf.config ~/.Arcitecta/mflux.cfg --csum-check --dest "$results_path" "$output_file_renamed"

echo "Finished processing sample: $sample"

#Log this job's resource usage stats

JOBID=$SLURM_JOB_ID

if [ ! -z $SLURM_ARRAY_JOB_ID ]; then

JOBID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

fi

my-job-stats -a -n -j $JOBID
