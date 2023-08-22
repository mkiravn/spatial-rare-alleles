# spatial-rare-alleles

Welcome to the project! 

You can run simulations by executing snakemake from the command line

`snakemake --cluster "sbatch -A pi-jnovembre -t 24:00:00 -p bigmem2 --mem=40000" --jobs 10`

Some of the larger simulations require running on bigmem2. The memory requirement for the simulation rule can be specified within it. Each simulation will end up in a folder in the `sims` directory with the name as the parameter combination.

After that, there are a couple more steps:
* Sample the SFS for each simulation run. You can do this in snakemake, or (as I've been doing more commonly) from the command line with `write_sfs_parallel.sh`. I've been running with 1000 samples and 100 iterations, but this can be changed within the script
* Now you need to compute the SFS into one file. You can do this with the command `Rscript scripts/compile.R inputdir` where `inputdir` is your directory of choice. This will create an output called `summary.df` with the compiled data.
