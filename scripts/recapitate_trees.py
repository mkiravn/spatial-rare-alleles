import pyslim as ps
import numpy as np
from scipy.stats import multivariate_normal as mvn
import tskit
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import msprime

np.random.seed(1234)
seed = np.random.randint(0,10000000)
input_file = "../1_W75_s-0.01_mu1e-10_K5_sigma0.2.trees"
num_to_sample = 1000

# Construct the file paths
# Get the directory and base filename from the input file path
file_dir = os.path.dirname(input_file)
base_filename = os.path.splitext(os.path.basename(input_file))[0]

# Create the output filename with the "sfs" directory included
num_samples_str = str(num_to_sample)
output_filename = os.path.join(file_dir, "sfs_log_niter", base_filename + '_' + num_samples_str + '.sfs.tsv')

# Load tree sequence
slim_ts = tskit.load(input_file)
slim_ts = ps.update(slim_ts)
locs = slim_ts.individual_locations  # locations of  all individuals

# Set the desired parameters
N_past = 1000  # Population size 2000 generations ago (2N)

r_ts = ps.sim_d(slim_ts,
            recombination_rate=1e-8,
            ancestral_Ne=N_past, random_seed=5)


# Initialize a list to store mutation times
mutation_times = []

# Iterate over mutations and collect their times
for mutation in slim_ts.mutations():
    mutation_times.append(10000 - mutation.metadata['mutation_list'][0]['slim_time'])

# Calculate summary statistics for the mutation times
mean_time = sum(mutation_times) / len(mutation_times)
min_time = min(mutation_times)
max_time = max(mutation_times)

print(f"Mean Mutation Age: {mean_time}")
print(f"Minimum Mutation Age: {min_time}")
print(f"Maximum Mutation Age: {max_time}")

# Create a histogram plot
plt.hist(mutation_times, bins=50, edgecolor='black')
plt.xlabel('Mutation Age')
plt.ylabel('Frequency')
plt.title('Mutation Age Histogram (s=0.01)')
plt.show()