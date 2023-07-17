import sys
sys.path.insert(0, '/home/mianniravn/.local/lib/python3.7/site-packages')
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
parser.add_argument('input_dir', help='Input directory containing tree sequence files')
parser.add_argument('--num_samples', type=int, default=100, help='Number of individuals to sample')
args = parser.parse_args()
print(args)

input_dir = args.input_dir
num_to_sample = args.num_samples
print(num_to_sample)
# Create output directory for SFS files
output_dir = os.path.join(input_dir, 'sfs',str(num_to_sample))
os.makedirs(output_dir, exist_ok=True)

print("Files to be analysed:")
print(os.listdir(input_dir))

# Process each file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('.trees'):
        # Construct the file paths
        input_path = os.path.join(input_dir, filename)
        num_samples_str =str(num_to_sample)
        output_filename = os.path.splitext(filename)[0] + '_' + num_samples_str + '.sfs.tsv'
        output_path = os.path.join(output_dir, output_filename)
        if os.path.exists(output_path):
               continue    
        print(input_path)
        # Load tree sequence
        slim_ts = tskit.load(input_path)
        slim_ts = ps.update(slim_ts)
        locs = slim_ts.individual_locations  # locations of  all individuals

        samps_at_end = ps.individuals_alive_at(slim_ts, 0)  # individuals alive at end
        
        def get_weights(lcs, mean=25, var=5):
            # gets sampling weights of individuals
            wts = np.zeros(len(lcs[:, 0]))
            sampkernel = mvn(mean=[mean, mean], cov=[[var, 0], [0, var]])

            for i in range(len(lcs[:, 0])):
                wts[i] = sampkernel.pdf(lcs[i, 0:2])

            wts /= np.sum(wts)  # normalize weights

            return wts


        def allele_counts(slim_ts, sample_sets=None):
            # returns allele counts in samples
            if sample_sets is None:
                sample_sets = [slim_ts.samples()]

            def f(x):
                return x

            return slim_ts.sample_count_stat(sample_sets=sample_sets, f=f, output_dim=len(sample_sets),
                                             span_normalise=False, windows='sites',
                                             polarised=True, mode='site', strict=False)


        # Gaussian samples with different variances
        vars_range = np.concatenate(([0.2], 0.2 * np.arange(0, 250, 10)))

        samples_dict = {}

        for var_value in vars_range:
            if var_value == 0:
                weights = np.ones(len(locs[:, 0])) / len(locs[:, 0])  # Uniform weights for random sample
            else:
                weights = get_weights(locs, 25, var_value)

            samples = np.random.choice(samps_at_end, size=num_to_sample, replace=False, p=weights)
            print(var_value)
            samples_dict[var_value / 0.2] = list(samples.astype(int))

        # Calculate SFS for each sample
        sfs_data = {}

        for var_value, samples in samples_dict.items():
            freqs = allele_counts(slim_ts, [samples])
            # convert the n x 1 array of floats to a vector of integers
            freqs = freqs.flatten().astype(int)
            mut_afs = np.zeros((len(samples) + 1, 1), dtype='int64')
            mut_afs[:, 0] = np.bincount(freqs, minlength=len(samples) + 1)
            sfs_data[var_value] = np.bincount(freqs, minlength=len(samples) + 1)

        # Create DataFrame from the SFS data
        sfs_df = pd.DataFrame(sfs_data)

        # Save SFS DataFrame to file
        sfs_df.to_csv(output_path, sep='\t', index_label='allele counts')


