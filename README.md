# UM-CoV-Seq

This repository holds pipelines, protocols, and web resources from the University of Montana to be used for the identification and characterization SARS-CoV-2 variants across the state of Montana, with an emphasis on serving rural and Tribal communities.

The main pipeline maps sequencing reads to the SARS-CoV-2 reference genome, then uses `GATK` and `iVar` to call variants and create consensus sequences. These consensus sequences are then assigned to lineages using `Pangolin` and `Nextclade`. The initial `iVar` pipeline was inspired by the pipeline from the Florida Coronavirus Genome Network (FL-CGNet) created at the University of Florida.

## Installation

### Clone repo

First, download this repo onto the computer where you will run these analysis:

    `git clone https://github.com/goodest-goodlab/UM-CoV-Seq.git`

Then, change working directories into this folder:

    `cd UM-CoV-Seq`

### Install software

This repo contains some of the software you need to run the pipeline (see Resources), but other programs need to be installed as well:

1. [snakemake](https://snakemake.readthedocs.io/en/stable/)- A pipeline management tool, which runs everything automatically.
2. [fastp](https://github.com/OpenGene/fastp)- For some processing of the raw reads (adaptor and quality trimming, merging of overlapping reads).
3. [bwa](http://bio-bwa.sourceforge.net)- For mapping reads.
4. [samtools](https://www.htslib.org)- For processing mapped SAM/BAM files.
5. [qualimap](http://qualimap.conesalab.org)- For QC of the mapped reads.
6. [MultiQC](https://multiqc.info)- For compiling QC reports across multiple samples.
7. [iVar](https://andersen-lab.github.io/ivar/html/index.html)- For variant calling and determining consensus sequences.
8. [GATK](https://gatk.broadinstitute.org/hc/en-us)- For variant calling.  
9. [bcftools](https://samtools.github.io/bcftools/)- For determining consensus sequences. 
10. [python](https://www.python.org)- For data processing.
11. [R](https://www.python.org)- For compiling results and making plots. 

You could install each program individually, but the easiest way to install them all is through the `conda` package manager from `Anaconda`. In the `envs` folder, we have provided the `UM-CoV-Seq.yml` pre-set environment to automatically install the necessary software. First, make sure Anaconda is installed on your system by following the instructions here: [Anaconda Download](https://www.anaconda.com/products/individual#download-section)

Then, start Anaconda:
    
    source /path/to/anaconda3/bin/activate

Load the environment:

    conda env create --file envs/UM-CoV-Seq.yml

And finally, activate the environment:

    conda activate UM-CoV-Seq

## Running the pipeline

### Specifying input data

Currently, `snakemake` is set up to work from a configuration file, stored in the `config_files` folder. You will need to create a new config file for each analysis batch. Config files are simple `.yaml` files that list parameters for the pipeline: which data to analyze, what reference genome to use, etc. Copy the provided template to make a new cofiguration file for your analysis.

### snakemake and SLURM

One advantage of `snakemake` is that it can handle job submission to the cluster for you. To do that, you need to create a `snakemake` "profile". A template is provided in the `profiles` folder. The profile is a single named folder, containing another file that must be named `config.yaml`. This file specifies how `snakemake` interacts with the SLURM job submission system used on Griz. Currently, you should only need to modify line 10 (to enter your email address to be notified if a job fails).

### Executing the pipeline. 

After you have installed all the necessary software, generated a config file for your analysis batch, and created a personal profile for submitting to the cluster, you are ready to run the pipeline. Make sure the `UM-CoV-Seq` environment is activated, and that the `UM-CoV-Seq` folder is your working directory. Then, execute:

```bash
snakemake -s pipe_reads_to_lineages.smk --configfile config_files/config_BATCH.yaml --profile profiles/YOUR_SNAKE_PROFILE_FOLDER/ --dryrun
```

The `--dryrun` flag means that `snakemake` will figure out what jobs it needs to do without actually running/submitting any of the jobs. It's good practice to do a dry run before you actually run the pipeline, to make sure it is doing what you think it should be doing. In particular, check that the number of total jobs matches your expectations: trimming and mapping steps need to happen for every sample, while most other steps should occur once. If all looks as it should, simplt remove the `--dryrun` flag to run the job for real:

```bash
snakemake -s pipe_reads_to_lineages.smk --configfile config_files/config_BATCH.yaml --profile profiles/YOUR_SNAKE_PROFILE_FOLDER/

```

And that should do it. As a final example, here's what the command looks like to analyze batch V4R1 using Tim's snakemake profile:

```bash
snakemake -s pipe_reads_to_lineages.smk --configfile config_files/config_V4R1.yaml --profile profiles/TJT/
```

In general, I would recommend running the pipeline within `screen` or `tmux`, as the pipeline can take a while (1-2 hours, depending on depth of sequencing and number of samples). The first run of the pipeline will be especially long, as some `R` packages need to be downloaded and installed. 

## Included resources

### SARS-CoV-2 reference genome

The reference genome is included in the folder `SARS-CoV-2-refseq`. It contains the full nucleotide assembly as well as the protein sequences and annotation information. 

This genome was downloaded from [NCBI Accession NC_045512](https://www.ncbi.nlm.nih.gov/genome/?term=NC_045512) on June 10, 2021.

More information can be found on [NCBI's SARS-CoV-2 Resources page](https://www.ncbi.nlm.nih.gov/sars-cov-2/)

### Masking of problematic sites

This repository also contains the [ProblematicSites_SARS-CoV2](https://github.com/W-L/ProblematicSites_SARS-CoV2) repository as a subtree, as there are several consistent errors in SARS-CoV-2 alignments that identify sites that should be removed from subsequent analyses (see: [Masking strategies for SARS-CoV-2 alignments](https://virological.org/t/masking-strategies-for-sars-cov-2-alignments/480)).

### Nextclade

Instead of putting it in the conda environment, We have also included the [NextClade binary](https://github.com/nextstrain/nextclade/releases) within the repo. 

# TODOs

1. Automatically generate GISAID output table. 
2. Decide on best practices for Pangolin versioning, and then implement. E.g., if we want to peg to a certain version, put that in the conda environment. Or, if we want to check for and update pangolin every time we reun pipeline, add in a rule or line of code to do that.
3. Update default SLURM profile: change location of SLURM log files.

# Ideas for the future

1. Run FastQC, or some other QC, immediately on the raw data? May be useful for diagnosing issues, though a lot of what we'd get from FastQC we probably already have from the Illumina software with the MiSeq.
2. More quality filtering during mapping? Could be more stringent in some spots. May be worth testing to see if it makes a difference.
3. Make use of the conda tools within `snakemake`? In theory, could set things up so that this whole pipeline has only 2 user-managed dependencies (`snakemake` and `conda`/`mamba`). `snakemkae` would handle the rest of the installing through `conda` environments for each rule. Might also improve things for running Nextclade.
4. Allow more flexibility in input file names? .fasta/.fa/.fastq/.fq, .gz or not? Right now, just works on files with `.fastq.gz` as their extension.
5. Improve the readgroup function? pull the first line of the FASTA and parse it to get a little more info, instead of the pretty simple stuff I'm doing currently.
6. Add a rule to plot the DAG?
7. Add a rule to output a .yml file of the conda environment used? 
