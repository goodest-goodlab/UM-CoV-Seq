# Tim's COVID analysis pipeline

I've made a basic pipeline that goes from raw `fastq.gz` reads and assigns lineages to each sample. This could be improved upon for the future (see some ideas in the TODOs), but seems to be orking fine for now. 

## Installing software for the pipeline

The pipeline uses a few tools:

1. [snakemake](https://snakemake.readthedocs.io/en/stable/)- A pipeline management tool, which basically runs everything automatically.
2. [fastp](https://github.com/OpenGene/fastp)- For some processing of the raw reads (adaptor and quality trimming, merging of overlapping reads).
3. [bwa](http://bio-bwa.sourceforge.net)- For mapping reads.
4. [samtools](https://www.htslib.org)- For processing mapped SAM/BAM files.
5. [qualimap](http://qualimap.conesalab.org)- For QC of the mapped reads.
6. [MultiQC](https://multiqc.info)- For compiling QC reports across multiple samples.
7. [iVar](https://andersen-lab.github.io/ivar/html/index.html)- For variant calling and determining consensus sequences.
8. [pangolin](https://cov-lineages.org/pangolin.html)- For assigning COVID sequences to lineages.

Right now, these tools can be installed using the [conda](https://docs.conda.io/en/latest/) package management tool, using the supplied environment description file.

You can create the environment with:

    conda env create --file envs/covid_pipeline_TJT.yml

And then activate the environment with:

    conda activate covid_pipeline_TJT

## Setting up the pipeline to work

### Specifying input raw reads

In the `pipe_reads_to_lineages.smk` file, you can modify the directory in which to search for raw data (hopefully self-explanatory in the file). The pipeline searches within the supplied directory (and any subdirectories) to find all `.fastq.gz` files that match the [Illumina file naming conventions](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm). Right now, the pipeline assumes that your data are:

1. Paired-end
2. Named according to the Illumina file conventions above.
3. All sample names are unique, and no samples have been sequenced multiple times (that is, there are only two files, one for R1 and one for R2, associated with each sample ID).

We could work on relaxing those assumptions in the future.

### Specifying the reference genome 

There's also a line in the `pipe_reads_to_lineages.smk` file for specifying the reference genome to map to and the `.GFF` file of annotations, which is again hopefully self-explanatory.

### Configuring snakemake

The `snakemake` script is set up to send most jobs to the SLURM job management system. The details of that are specified in the `config.yaml` file in the `snake_profile_TJT` folder.  You'll need to make your own `snakemake` profile folder, copy the `config.yaml` file into it, and then make a few changes:

1. I have it set up to output SLURM logs to a folder within my copy of the `UM-CoV-Seq` repo. You'll want to substitute in your own folder. NOTE- stupidly, the folder must exist already, SLURM can't/won't make the directory if it doesn't already exist and the log files will disappear into the ether. So ensure that the directory exists. 
2. I have it set up to email you if a job fails. Put in your email address, so it doesn't email me!

## Running the pipeline

If you've set up the pipeline as listed above, invoking the pipeline to run should be very easy. Make sure the `covid_pipeline_TJT` conda environment is activated. `conda` can get weird about which version of Python it is using, run `which python` and make sure that it is using the Python install within the `covid_pipeline_TJT/env/` folder. If not, you may need to deactivate and reactivate the environment. 

Once that's all set up, just navigate to the `UM-CoV-Seq` directory that you're working from, and execute (using your own snake profile folder):

```bash
snakemake -s pipe_reads_to_bams.smk --profile YOUR_SNAKE_PROFILE_FOLDER/ --dryrun
```

The `--dryrun` flag means that `snakemake` will figure out what jobs it needs to do without actually running/submitting any of the jobs. It's good practice to do a dry run before you actually run the pipeline, to make sure it is doing what you think it should be doing. If it is, just remove the `--dryrun` flag:

```bash
snakemake -s pipe_reads_to_bams.smk --profile YOUR_SNAKE_PROFILE_FOLDER/ 
```

And that should do it. As a final example, here's what the command looks like using my (Tim's) profile:

```bash
snakemake -s pipe_reads_to_bams.smk --profile snake_profile_TJT/ 
```

I would recommend running the pipeline within `screen` or `tmux` or whatever, as the whole thing could take a while. You can always use `nohup ... &` as well, I guess, but I like something like `tmux` better. 

# TODOs

1. Figure out whether to add the -aa option for the pileup for consensus sequences for iVar. 

# IDEAS

1. Run FastQC, or some other QC, immediately on the raw data?
2. more quality filtering during mapping?
3. Make use of the conda tools within `snakemake`? In theory, could set things up so that this whole pipeline has only 2 user-managed dependencies (`snakemake` and `conda` or `mamba`). `snakemkae` would handle the rest of the installing through `conda` environments. 
4. Allow more flexibility in input file names? .fasta/.fa/.fastq/.fq, .gz or not? Right now, just works on files with `.fastq.gz` as their extension.
5. Mark all intermediate files as temp in snakemake, to save disk space? Not that useful now, maybe once it is finalized?
6. Improve Nextclade? Edit default for QC, use our own .gff (which didn't work initially), etc.
7. Improve the readgroup function? pull the first line of the FASTA and parse it to get a little more info, instead of the pretty simple stuff I'm doing currently.
8. Add a rule to plot the DAG?
