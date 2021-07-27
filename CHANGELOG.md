2021-07-19 Timothy Thurman <timothy.j.thurman@gmail.com>

    * Added this CHANGELOG.

    * Updated profile folder with template `config.yaml` file, changed naming and location of slurm logs. 

    * Introduce config files, so snakemake pipeline no longer needs to be edited. 

    * Some minor updates to snakefile requested by linter. (1) all thread/CPU usge specified by slurm variables instead of snakemake resources (2) added log files for all rules.

    * Pipeline now outputs elapsed time and preserves copy of the snakefile used.

    * Add a check that the data folder exists, for sensible snakemake errors if it doesn't.

    * Added my old .yml files to the envs folder. Tested a new, from-scratch conda environment (UM-CoV-Seq.yml).

    * iVar variants and iVar consensus sequence now use same mpileup command. 

    * Merge READMEs and snakemake pipelines, so we don't have TJT and GWCT copies anymore. 


2021-07-18 Gregg Thomas <greggwct@gmail.com>

    * Changes currently listed by TJT, Gregg can add more in later.

    * A major update. Gregg added his GATK/bcftools variant calling and consensus sequence steps into the snakemake pipline, as pipe_read_to_lineage_gwct.smk. Includes new step to mask problematic sites out of the pileup before variant calling or consensus seq making with iVar, and also masks .vcfs for consensus sequences from bcftools. 

    * In new pipeline, now only need to specify the batch name provided by David for a given sequencing batch, finds the corresponding raw data folder and uses the batch name for the output data folders. 

    * New scripts, in `lib`, for masking of vcfs and pileups, counting bases in consensus seq files (removes some dependencies), and formatting sample IDs for GISAID uploading.

    * compile_results.R: Updated to remove some dependencies, self-manage the remaining ones, and do some plotting as well as compilation of results.

    * Removes renv stuff, no longer needed.    

2021-06-21  Timothy Thurman  <timothy.j.thurman@gmail.com>

    * pipe_read_to_lineage.smk: Added step to compile QC, variant number, pangolin, and Nextclade results into a summary. 

    * Setup renv to manage the R dependencies needed to compile the final results. 

2021-06-21  Timothy Thurman  <timothy.j.thurman@gmail.com>

    * CHANGELOG.md - Made this changelog to track changes

    * Added `nextclade-resources`, with default tree, qc, reference, and gff files for running the same (default) Nextclade analysis you get if drag-and-drop a file online. 

    * pipe_read_to_lineage.smk: All output files now batched into user-specified directory, to isolate acros sequencing runs or when using different parameters. Added a preliminary Nextclade output. Added a masking step, but need to run with/without to see whether it is really doing anything. 