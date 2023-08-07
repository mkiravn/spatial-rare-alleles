### Welcome to the spatial rare alleles simulation pipeline!
### This pipeline runs spatial simulations in SLiM
### and then samples individuals with a spatial kernel
### to obtain an SFS


replicates = 100 # number of simulation replicates
s_coeffs = [-1e-3, -1e-2, -1e-1]
s_coeffs_slim =  [s*2 for s in s_coeffs] # selection coefficients for slim
mus = [1e-10] # mutation rates
Ks = [5] # densities
Ws = [75]
sigmas = [0.2]

rule all:
    input:
        expand("sfs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n1000_niter1000.sfs",
               rep=range(replicates),
               s=s_coeffs,
               mu=mus,
               K=Ks,
               W=Ws,
               sigma=sigmas)


# Define the first step: run_simulations
rule run_simulations:
    output:
        "sims/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees.log"
    params: s_val= lambda wildcards: str(wildcards.s) 
    resources: mem_mb=40000
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
        "sims/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.trees"
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.sfs.log"
    output:
        "sfs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}_n{n}_niter{niter}.sfs"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {wildcards.n} \
            --output {output} \
            --niter {wildcards.niter}
        """
