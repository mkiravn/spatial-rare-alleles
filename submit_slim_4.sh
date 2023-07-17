#!/bin/bash

#SBATCH --job-name=space_simu
#SBATCH --output=logs/space_simu_%A_%a.out
#SBATCH --error=logs/space_simu_%A_%a.err
#SBATCH --array=0-3%5 # changed
#SBATCH --time=03:00:00  # 4 days
#SBATCH --partition=bigmem2
#SBATCH --account=pi-jnovembre
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  # Use 1 CPU core per task
#SBATCH --mem-per-cpu=50G

arr_s=(-2e-3 -2e-2 -2e-1)
arr_mu=(1e-10)
replicates=100
W=100  # Fixed value for W

timestamp=$(date +"%Y%m%d%H%M")  # Get the current timestamp
output_file="slim_times_sfs_$timestamp.txt"  # Append timestamp to the output file name
mkdir -p sims/$timestamp

# Define the index variables for array indexing
s_index=$(($SLURM_ARRAY_TASK_ID % ${#arr_s[@]}))
mu_index=$(($SLURM_ARRAY_TASK_ID / ${#arr_s[@]} % ${#arr_mu[@]}))

s=${arr_s[$s_index]}
mu=${arr_mu[$mu_index]}
replicate=$((SLURM_ARRAY_TASK_ID % replicates + 1))

echo "Combination: s=$s, mu=$mu, W=$W, replicate=$replicate, timestamp=$timestamp"

module load SLiM/3.0 

# Function to run a single simulation
run_simulation() {
    start_time=$(date +%s.%N)
    echo "Command will be slim -seed $((1234 + $replicate)) -time -mem -long -define sigma=0.2 -define mu=$mu -define s=$s -define W=$W -define id=$replicate -define stamp=$timestamp"
    slim -seed $((1234 + $replicate)) -time -mem -long -define sigma=0.2 -define mu=$mu -define s=$s -define W=$W -define id=$replicate -define stamp=$timestamp recipe_space.slim &> "logs/slim_output.txt"
    end_time=$(date +%s.%N)
    duration=$(echo "$end_time - $start_time" | bc)
    echo "Time taken: $duration seconds"
    echo "$s $mu $replicate $W $duration" >> "$output_file"
}

# Run the simulation
run_simulation &
wait

