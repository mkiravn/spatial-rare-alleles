### Welcome to the spatial rare alleles simulation pipeline!
### This pipeline runs spatial simulations in SLiM
### and then samples individuals with a spatial kernel
### to obtain an SFS


replicates = 2 # number of simulation replicates
s_coeffs = [-1e-3, -1e-2, -1e-1]
s_coeffs_slim =  [s*2 for s in s_coeffs] # selection coefficients for slim
mus = [1e-10] # mutation rates
Ks = [5] # densities
Ws = [75]
sigmas = [0.2]

rule all:
    input:
        expand("sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_N{past_Ne}_t{time_bottle}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_recapitated_N{past_Ne}_t{time_bottle}.trees",
               rep=range(replicates),
               s=s_coeffs,
               mu=mus,
               K=Ks,
               W=Ws,
               sigma=sigmas,
		past_Ne=1000,time_bottle=2000)

# Define the first step: run_simulations
rule run_simulations:
    output:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees.log"
    params: s_val= lambda wildcards: str(wildcards.s) 
    resources: mem_mb=60000
    shell:
        """
        module load SLiM/3.0
        slim -define K={wildcards.K} \
            -define mu={wildcards.mu} \
            -define s={params.s_val} \
            -define W={wildcards.W} \
            -define "outpath='{output}'" \
            -define sigma={wildcards.sigma} \
            recipe_space_sm.slim
         """

# Second step: sfs sampling
rule sample_sfs:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.sfs.log"
    output:
        "sfs/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.sfs"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {wildcards.n} \
            --output {output} \
            --niter {wildcards.niter}
        """

# Variant sampling
rule sample_variants:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.variants.log"
    output:
        "sfs/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.variants"
    shell:
        """
        python scripts/parallel_sfs_maggie.py {input} \
            --num_samples {wildcards.n} \
            --output {output} \
            --niter {wildcards.niter}
        """


# checking spatial positions
rule sample_locations:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}.sds.log"
    output:
        "sds/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}.sds"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {wildcards.n} \
            --output {output}
        """

# cut and recapitate trees
rule cut_recapitate:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_recapitated_N{past_Ne}_t{time_bottle}.trees.log"
    output:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_N{past_Ne}_t{time_bottle}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_recapitated_N{past_Ne}_t{time_bottle}.trees"
    shell:
        """
        python scripts/cut_recapitate_snakemake.py {input} \
            --past_Ne {wildcards.past_Ne} \
            --time_bottle {wildcards.time_bottle} \
            --output {output}
        """

# Second step: sfs sampling
rule sample_sfs_recapitated:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}_recapitated_N{past_Ne}_t{time_bottle}.sfs.log"
    output:
        "sfs/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_N{past_Ne}_t{time_bottle}_n{n}_niter{niter}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_recapitated_N{past_Ne}_t{time_bottle}.sfs"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {wildcards.n} \
            --output {output} \
            --niter {wildcards.niter}
        """


# Variant sampling
rule sample_variants_recapitated:
    input:
        "sims/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}_recapitated_N{past_Ne}_t{time_bottle}.variants.log"
    output:
        "sfs/W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_N{past_Ne}_t{time_bottle}_n{n}_niter{niter}/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_recapitated_N{past_Ne}_t{time_bottle}.variants"
    shell:
        """
        python scripts/parallel_sfs_maggie.py {input} \
            --num_samples {wildcards.n} \
            --output {output} \
            --niter {wildcards.niter}
        """
