Code and data supplement to the manuscript titled "Rapidly changing speciation and extinction rates can be inferred in spite of non-identifiability", by Kopperud, Magee & Höhna

## Organization

Scripts are run from the top level. For example R scripts are run as: `Rscript scripts/hypothetical_models_frommu.R`, and bash scripts are run as `./scripts/create_SLURM_JOBS.sh` in the terminal.

The sub-directories are:

* `trees`: the tree files
* `scripts`: the scripts that we used, in R, bash and RevBayes
* `figures`: the figures produced by the scripts are stored here
* `jobs`: the SLURM scripts
* `logs`: the RevBayes screen output
* `output`: the RevBayes log files, i.e. the posterior samples

## 01: RevBayes episodic birth-death rate estimation

We generate a set of SLURM scripts that we run on the computer cluster using `bash scripts/create_SLURM_JOBS.sh`. We submit these jobs to the SLURM service using the command `sbatch`. The result of the episodic birth-death models (the trace files) are placed in the directory `output`. We did not commit these since they are on the order of 25 GB.

The revbayes script `analysis.Rev` takes three arguments: the tree file, the number of repeated MCMC runs, and the seed for reproducibility.

The `analysis.Rev` script takes care of the MCMC details, but calls the `HSMRBDP.Rev` script which lays out how the episodic birth-death model is set up.

## 02: Convergence assessment

The `convergence_plots.R` script performs the Kolmogorov-Smirnov test and produces the plots that are shown in the supplementary material of the manuscript.

## 03: Computing the posterior median and credible interval

We use the script under `scripts/plot_HSMRF.R` to 1) process the trace files and compute quantiles for the rates for each interval, 2) save those summaries to a file, and 3) to plot the resulting models so that we can inspect the results visually.

## 04: Analysis of congruence classes

These scripts analyze the congruence class and produce the figures presented in the ms:
* `empirical_congruence_classes.R`: the empirical datasets
* `hypothetical_models_frommu.R`: the hypothetical scenarios where we proposed alternative extinction rates
* `hypothetical_models_fromlambda.R`: the hypothetical scenarios where we proposed alternative speciation rates

