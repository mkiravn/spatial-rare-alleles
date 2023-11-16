#### A script which takes a recapitated tree and titrates in neutral mutations, and finally samples the SFS

# import modules
import sys
sys.path.append('scripts')
import pyslim as ps
import numpy as np
import tskit
import argparse
import os
import msprime
from sfs_function import *

np.random.seed(1234)
gen_time = 5.6 # from simulations


input_file = "../sims/1_W75_s-0.01_mu1e-10_K5_sigma0.2.trees"
past_Ne = 1000
time_bottle = 2000


# Load tree sequence
slim_ts = tskit.load(input_file)
slim_ts = ps.update(slim_ts)

# compute mutation times
tables = slim_ts.dump_tables()
tables.compute_mutation_times()
tables.sort()
slim_ts = tables.tree_sequence()



# specify ancestry:
N_past = past_Ne


# mutate
rates = [-10]
times=[500,1000,5000]
sizes = [500,1000,5000,10000]
sfs_df=pd.DataFrame()
variants_df=pd.DataFrame()

for rate in rates:
    for time_bottle in times:
        for N_past in sizes:
            print(f"Sampling for rate {rate}, size {N_past} and time {time_bottle}")
            # we decapitate
            new_ts = slim_ts.decapitate(2000 * gen_time)
            # construct demography
            demography = msprime.Demography.from_tree_sequence(slim_ts,
                                                               initial_size=slim_ts.num_individuals)

            # adjust population metadata for demographic model
            i = 0
            for pop in demography.populations:
                # must set their effective population sizes
                pop.initial_size = slim_ts.num_individuals * i * gen_time
                pop.initially_active = i
                i += 1

            # Population size change from 2000 generations ago to the beginning of the simulation
            demography.add_population_parameters_change(time=time_bottle * gen_time,
                                                        # in past. Scaling needed im order to match SLiM time
                                                        initial_size=N_past * gen_time,  #
                                                        population="pop_1")

            # recapitate tree sequence
            new_ts = msprime.sim_ancestry(initial_state=new_ts,
                                          recombination_rate=1e-8 / gen_time,
                                          demography=demography)  # population halves in size

            print(f"rate: {1 * 10**(rate)}")
            mut_ts = msprime.sim_mutations(new_ts,
                                   rate = 1 * 10**(rate) / gen_time, # fewer mutations per tick since 1 tick ~= 1/5 generations
                                   model=msprime.SLiMMutationModel(type=0, slim_generation = 10000), # time ensures all the mutation times match up
                                   keep=True,
                                   random_seed=7)
            propneut=(np.shape(mut_ts.tables.mutations.metadata)[0]-np.shape(slim_ts.tables.mutations.metadata)[0])/(np.shape(mut_ts.tables.mutations.metadata)[0])
            print(f"Proportion of neutral mutations = {propneut}")
            sfs_i = sfs(mut_ts,n_bins=20,niter=50,num_to_sample=10000,output=f"{input_file}_{rate}.sfs",sd_range_disp=np.array([10,100]))
            sfs_i = sfs_i.assign(rate=rate,propneut=propneut,N_past=N_past,time_bottle=time_bottle)
            sfs_df = pd.concat([sfs_df,sfs_i],axis=0)

            variants_i = biplot(mut_ts,niter=50,num_to_sample=10000,sd_range_disp=np.array([10,100]),subsample_loci=False)
            variants_i = variants_i.assign(rate=rate,propneut=propneut,N_past=N_past,time_bottle=time_bottle)
            variants_df = pd.concat([variants_df, variants_i], axis=0)
sfs_df.to_csv(f"{input_file}_allpars.sfs", sep='\t', index=False, na_rep="NA")
variants_df.to_csv(f"{input_file}_allpars.variants", sep='\t', index=False, na_rep="NA")