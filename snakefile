


replicates = 100 # number of simulation replicates
s_coeffs = [-1e-3, -1e-2, -1e-1]
s_coeffs_slim =  [s*2 for s in s_coeffs] # selection coefficients for slim
mus = [1e-10] # mutation rates
Ks = [5] # densities
Ws = [75]
sigmas = [0.2]

rule all:
    input:
        expand("sfs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.sfs",
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
    params:
        n = 1000,
        niter=1000
    log: "logs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.sfs.log"
    output:
        "sfs/{rep}_W{W}_s{s}_mu{mu}_K{K}_sigma{sigma}.sfs"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {params.n} \
            --output {output} \
            --niter {params.niter}
        """
