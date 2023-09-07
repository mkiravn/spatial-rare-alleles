#!/bin/bash

#SBATCH --job-name=space_sfs
#SBATCH --output=logs/space_sfs_%A_%a.out
#SBATCH --error=logs/space_sfs_%A_%a.err
#SBATCH --array=0-99%3
#SBATCH --time=00:45:00
#SBATCH --partition=broadwl
#SBATCH --account=pi-jnovembre
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  # Use 1 CPU core per task
#SBATCH --mem=5G

conda activate
module unload python
module load python/cpython-3.8.5
#module load anaconda-2021.05
#source activate


# Add lines here to run your Python script on each task
file_dir="sims"
file_list=($file_dir/*.trees)
filename=${file_list[$SLURM_ARRAY_TASK_ID-1]}
output_dir="sfs_new/"
output_file="${filename}_n1000_niter100.variants.tsv"

# Function to run a single simulation
run_analysis() {
    start_time=$(date +%s.%N)
    # Check if the output file already exists
    if [ -f "$output_dir$output_file" ]; then
      echo "Output file $output_dir$output_file already exists. Skipping processing."
    else
      echo "Processing input file $filename"
      # Run the Python script
      python scripts/parallel_sfs_maggie.py --num_samples 1000 --niter 100 --output $output_dir$output_file $filename
    fi
    echo "Time taken: $duration seconds"
}

# Run the simulation
run_analysis &
wait

