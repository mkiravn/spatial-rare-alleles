import sys
#sys.path.insert(0, '/home/mianniravn/.local/lib/python3.7/site-packages')
from matplotlib import pyplot as plt
import pyslim as ps
import numpy as np
from scipy.stats import multivariate_normal as mvn
import tskit
import pandas as pd
import argparse
import os

np.random.seed(1234)

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Calculate and save SFS data.')
parser.add_argument('input_file', help='File to be processed')
parser.add_argument('--num_samples', type=int, default=100, help='Number of individuals to sample')
parser.add_argument('--output', help='Output file path')
args = parser.parse_args()
print(args)

input_file = args.input_file
num_to_sample = args.num_samples
output = args.output
print(num_to_sample)

num_bins = 50 # number of bins for bincount


# Construct the file paths
# Get the directory and base filename from the input file path
file_dir = os.path.dirname(input_file)
base_filename = os.path.splitext(os.path.basename(input_file))[0]

# Create the output filename with the "sds" directory included
num_samples_str =str(num_to_sample)

# Load tree sequence
slim_ts = tskit.load(input_file)
slim_ts = ps.update(slim_ts)
locs = slim_ts.individual_locations  # locations of  all individuals

samps_at_end = ps.individuals_alive_at(slim_ts, 0)  # individuals alive at end

# genotype matrix - actually a haplotype matrix
H = slim_ts.genotype_matrix(isolated_as_missing=False)

# Assuming G is a numpy array of shape (L, 2N)
L, N = H.shape[0], H.shape[1] // 2

# Reshape G to (L, N, 2) to separate the two genomes
G_reshape = np.reshape(H, (L, N, 2))

# Sum the two genomes in each individual to get a genotype matrix
G = np.sum(G_reshape, axis=2)
G[G<1]=0 # fixes an occasional mysterious -1 which appears




def get_weights(lcs, mean=25, var=5):
    # gets sampling weights of individuals
    wts = np.zeros(len(lcs[:, 0]))
    sampkernel = mvn(mean=[mean, mean], cov=[[var, 0], [0, var]])

    for i in range(len(lcs[:, 0])):
        wts[i] = sampkernel.pdf(lcs[i, 0:2])

    wts /= np.sum(wts)  # normalize weights

    return wts


# Gaussian samples with different sds
sd_range = np.concatenate((np.arange(0,2,0.1), np.arange(5, 105, 5)))


samples_dict = {}
sd_list =[]
for sd_value in sd_range:
    # var_value = (sd_value * np.sqrt(0.06))**2 # re-scaling for sigma=0.2 and seff=sqrt(0.06) 
    var_value = sd_value **2
    if var_value == 0:
        weights = np.ones(len(locs[:, 0])) / len(locs[:, 0])  # Uniform weights for random sample
    else:
        weights = get_weights(locs, 25, var_value)

    samples = np.random.choice(samps_at_end, size=num_to_sample, replace=True, p=weights)

    samples_dict[sd_value] = list(samples.astype(int))

print(f"Computed weights and chose samples for {input_file}")

# Calculate SFS for each sample
sfs_data = {}

bin_edges = np.geomspace(1,num_to_sample,num=num_bins)

for sd_value, samples in samples_dict.items():
    empi_sd_x = np.std(locs[samples,0])
    print(f"Sampling sd: {sd_value} Sd of locations: {empi_sd_x}")
    sd_list.append(empi_sd_x)

# Create DataFrame from the SFS data
sfs_data['allele counts'] = bin_edges[1:len(bin_edges)]
print(len(sfs_data['allele counts']))
sfs_df = pd.DataFrame(sfs_data)

# Save SFS DataFrame to file
print("Saving file to")
print(output)

sds = {"Sampling":sd_range,"Realised":sd_list}
sds_df = pd.DataFrame(sds)
sds_df.to_csv(output, sep='\t', index=False,na_rep= "NA")




