#!/bin/bash

#SBATCH --job-name=get_sfs
#SBATCH --output=logs/get_sfs_%A_%a.out
#SBATCH --error=logs/get_sfs_%A_%a.err
#SBATCH --time=24:00:00  # 1 day
#SBATCH --partition=jnovembre
#SBATCH --account=pi-jnovembre
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1  # Use 1 CPU core per task
#SBATCH --mem-per-cpu=5G

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"

# Load required modules
module load python
conda activate


# Add lines here to run your Python script on each task
file_dir="sims/20230701/"
file_list=($file_dir*.trees)
filename=${file_list[$SLURM_ARRAY_TASK_ID-1]}


output_dir="sims/20230701/sfs/1000/"
output_file="${filename}_1000.sfs.tsv"

# Check if the output file already exists
if [ -f "$output_dir$output_file" ]; then
  echo "Output file $output_dir$output_file already exists. Skipping processing."
else
  echo "Processing input file $filename"
  # Run the Python script
  python scripts/parallel_sfs.py $filename --num_samples 1000
fi

