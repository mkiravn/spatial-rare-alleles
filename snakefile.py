


replicates = 100 # number of simulation replicates
s_coeffs = [-1e-3 -1e-2 -1e-1]
s_coeffs_slim =  [s/2 for s in s_coeffs] # selection coefficients for slim
mus = [1e-10] # mutation rates
Ks = [5] # densities

rule all:
    input:
        expand("sfs/{rep}_W{w}_s{s}_mu{mu}_K{K}.sfs",
               rep=range(replicates),
               s=s_coeffs_slim,
               mu=mus,
               K=Ks)


# Define the first step: run_simulations
rule run_simulations:
    output:
        "sims/{rep}_W{w}_s{s}_mu{mu}_K{K}.trees"
    log: "logs/{rep}_W{w}_s{s}_mu{mu}_K{K}.trees.log"
    shell:
        """
        slim -define K={wildcards.K} \
            -define mu={wildcards.mu} \
            -define s={wildcards.s} \
            -define W={wildcards.W} \
            -define outpath={output} \
            recipe_space_sm.slim
         """

# Second step: sfs sampling
rule sample_sfs:
    input:
        "sims/{rep}_W{w}_s{s}_mu{mu}_K{K}.trees"
    params:
        n = 1000,
        niter=1000
    log: "logs/{rep}_W{w}_s{s}_mu{mu}_K{K}.sfs.log"
    output:
        "sfs/{rep}_W{w}_s{s}_mu{mu}_K{K}.sfs"
    shell:
        """
        python scripts/parallel_sfs_niter.py {input} \
            --num_samples {params.n} \
            --output {output} \
            --niter {params.niter}
        """
