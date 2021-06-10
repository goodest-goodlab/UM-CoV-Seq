# UM-CoV-Seq

This repository will hold pipelines, protocols, and web resources from the University of Montana to be used for the identification and characterization SARS-CoV-2 variants across the state of Montana, with an emphasis on serving rural and Tribal communities.

The initial pipeline comes from the Florida Coronavirus Genome Network (FL-CGNet) created at the University of Florida.

This repository also contains the [ProblematicSites_SARS-CoV2](https://github.com/W-L/ProblematicSites_SARS-CoV2) repository as a submodule, as there are several consistent errors in SARS-CoV-2 alignments that identify sites that should be removed from subsequent analyses (see: [Masking strategies for SARS-CoV-2 alignments](https://virological.org/t/masking-strategies-for-sars-cov-2-alignments/480)).

## Installing pipeline as a conda environment

This pipeline uses genomics software to go from raw SARS-CoV-2 sequence data to lineage-assigned phylogenies of samples. This software includes:

1. [BWA](http://bio-bwa.sourceforge.net/) for read mapping.
2. [iVar](https://andersen-lab.github.io/ivar/html/) for variant calling.
3. The aforementioned [ProblematicSites_SARS-CoV2](https://github.com/W-L/ProblematicSites_SARS-CoV2) and associated scripts for site masking.
4. [IQ-Tree](http://www.iqtree.org/) for phylogeny inference.
5. [Pangolin](https://github.com/cov-lineages/pangolin) for lineage assignment.

You will need each of these programs installed on your system to run the pipeline. You can install them each individually, but they can also easily be installed from [bioconda](https://anaconda.org/bioconda), and we have provided the `UM-CoV-Seq.yml` pre-set environment to automatically install them.

First, make sure Anaconda is installed on your system by following the instructions here: [Anaconda Download](https://www.anaconda.com/products/individual#download-section)

Then, start Anaconda:
    
    source /path/to/anaconda3/bin/activate

Load the environment:

    conda env create --file UM-CoV-Seq.yml

And finally, activate the environment:

    conda activate UM-CoV-Seq

### More to come!

