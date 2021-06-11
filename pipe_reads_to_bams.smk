############################
##   EDIT BEFORE RUNNING  ##
############################

# List the directory containing the raw reads you 
# want to process:
RAW_DATA_FOLDER="/home/tt164677e/tim_beegfs/raw_read_data/SARS-CoV-2-06102021-Illumina"

# List the reference genome you want to map to:
REF="/home/tt164677e/tim_beegfs/UM-CoV-Seq/SARS-CoV-2-refseq/GCF_009858895.2_ASM985889v3_genomic.fna"

#################################
## DO NOT EDIT BELOW THIS LINE ##
##  WHEN RUNNING THE PIPELINE  ##
#################################


###############
##   SETUP   ##
############### 

# Python env setup
import os
import re

# Pull all sample files from the raw data folder
# get filenames of the individual fastas
fasta_fullpath = []
samples = []
for root, dirs, files in os.walk(RAW_DATA_FOLDER):
    for name in files:
        if re.search("_S\d+_L\d+_R[12]_001.fastq.gz$", name):
            samples.append(re.sub("_S\d+_L\d+_R[12]_001.fastq.gz$", "", name))
            fasta_fullpath.append(os.path.join(root, name))

# Get unique sample IDs
samples = list(set(samples))
samples.sort()

# Get filename for the genome index
index_path= REF + ".amb"


######################
## HELPER FUNCTIONS ##
######################

# At the very start, need to match each
# sample ID back to its fastq files.
# This assumes that each sample only has one set of sequence files in a given folder. 
def get_R1_for_sample(wildcards):
    outfile = "file_not_found.txt"
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("R1", filename):
                outfile = filename
    return outfile

def get_R2_for_sample(wildcards):
    outfile = "file_not_found.txt"
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("R2", filename):
                outfile = filename
    return outfile



# A function to create a readgroup for BWA from the sample, lane, and run info
def make_RG(wildcards):
    for filename in fasta_fullpath:
        if re.search(wildcards.sample, filename):
            if re.search("R1", filename):
                basename = os.path.basename(filename)
    # Extract sample ID, lane, and run (seq1,seq2, seq3) from input
    sample_ID = basename.split("_")[0]
    sample_num = basename.split("_")[1]
    lane = basename.split("_")[2]
    # Assemble the RG header. Fields:
    # ID: Individual sample ID plus sample number
    # LB: library, sample + "lib1"
    # PL: platform, ILLUMINA
    # SM: sample, sample
    # PU: platform unit, run + lane + sample
    rg_out = "@RG\\tID:" + sample_ID + sample_num + "\\tLB:" + sample_ID + "\\tPL:ILLUMINA" + "\\tSM:" + sample_ID
    return rg_out

####################
## PIPELINE START ##
####################
localrules: all

rule all:
    input:
        "processed_reads/QC/multiqc/multiqc_report_trimming.html", # QC report on trimming
        "processed_reads/QC/multiqc/multiqc_report_dedup_Picard.html", # QC report on dedup process
        "processed_reads/QC/multiqc/multiqc_report_raw_bams.html", # QC report on raw bams
        "processed_reads/QC/multiqc/multiqc_report_dedup_bams.html" # QC report on final bams


## trim_raw_reads : remove adaptors and low-quality bases
# Uses fastp. Options:
# -m Merge mode: merge overlapping read pairs which overlap
# -c Correct mismatched bases in the merged region
# Using default parameters for merge and correction
# --detect_adapter_for_pe Auto-detects possible adaptor sequences for removal
# --cut_front Do a sliding window analysis from the 5' end, cut read when quality falls below thresh
# --cut_front_window_size Window size of 5 for sliding window 
# --cut_front_mean_quality Mean quality of 20 for sliding window
# -l 25 Minimum length of 25 BP for the read
# -w Number of cores
# -h, -j Name of report HTML and JSON  output reports
rule trim_and_raw_reads:
    input:
        raw_r1=get_R1_for_sample,
        raw_r2=get_R2_for_sample
    output:
        trim_merged="processed_reads/trimmed/{sample}.merged.fq.gz",
        trim_r1_pair="processed_reads/trimmed/{sample}.nomerge.pair.R1.fq.qz",
        trim_r2_pair="processed_reads/trimmed/{sample}.nomerge.pair.R2.fq.qz",
        trim_r1_nopair="processed_reads/trimmed/{sample}.nopair.R1.fq.qz",
        trim_r2_nopair="processed_reads/trimmed/{sample}.nopair.R2.fq.qz",
        rep_html="logs/fastp/{sample}_trim_fastp.html",
        rep_json="logs/fastp/{sample}_trim_fastp.json"
    resources:
        cpus = 4
    # conda:
    #     "envs/UM_COV_TJT.yml"
    log:
        log="logs/fastp/{sample}_trim_log.txt"
    shell:
        """
        fastp -i {input.raw_r1} -I {input.raw_r2} -m --merged_out {output.trim_merged} --out1 {output.trim_r1_pair} --out2 {output.trim_r2_pair} --unpaired1 {output.trim_r1_nopair} --unpaired2 {output.trim_r2_nopair} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.rep_json} -h {output.rep_html} -w {resources.cpus} 2> {log.log}
        """

## multiqc_trim_reports: collate fastp trimming reports
rule multiqc_trim_reports:
    input:
        expand("logs/fastp/{sample}_trim_fastp.json", sample = samples)
    output:
        "processed_reads/QC/multiqc/multiqc_report_trimming.html"
    shell:
        """
        multiqc logs/fastp -o processed_reads/QC/multiqc/ -n multiqc_report_trimming.html
        """

## index_ref: index genome for BWA
# Index the reference genome, if it isn't already
rule index_ref:
    input:
        REF
    output:
        multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")    
    shell:
        """
        bwa index {input}
        """

## map_merged_reads: map trimmed, merged reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_merged_reads:
    input:
        reads="processed_reads/trimmed/{sample}.merged.fq.gz",
        genome=REF,
        genome_index=index_path
    output:
        "processed_reads/mapped/{sample}.merged.sorted.bam"
    params:
        basename="processed_reads/mapped/{sample}",
        read_group=make_RG
    log:
        "logs/mapping/{sample}_merged.log"
    resources:
        cpus=8
    shell:
        """
        # Run bwa mem, pipe to samtools view to convert to bam, save as a tmp.bam
        bwa mem -M -t {resources.cpus} -R '{params.read_group}' {input.genome} {input.reads} 2> {log} | samtools view -b - > {params.basename}.merged.unsorted.tmp.bam

        # then sort
        samtools sort {params.basename}.merged.unsorted.tmp.bam -o {output}

        # Then remove the intermediate, unsorted .bam
        rm {params.basename}.merged.unsorted.tmp.bam
        """

# # map_unmerged_pairs: map trimmed, not merged, paired reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_unmerged_pairs:
    input:
        reads_forward="processed_reads/trimmed/{sample}.nomerge.pair.R1.fq.qz",
        reads_reverse="processed_reads/trimmed/{sample}.nomerge.pair.R2.fq.qz",
        genome=REF
    output:
        "processed_reads/mapped/{sample}.nomerge.paired.sorted.bam"
    params:
        basename="processed_reads/mapped/{sample}",
        read_group=make_RG
    log:
        "logs/mapping/{sample}_nomerge_paired.log"
    resources:
        cpus=24
    shell:
        """
        # Run bwa mem, pipe to samtools view to convert to bam, save as a tmp.bam
        bwa mem -M -t {resources.cpus} -R '{params.read_group}' {input.genome} {input.reads_forward} {input.reads_reverse} 2> {log} | samtools view -b - > {params.basename}_nomerge.paired.unsorted.tmp.bam

        # then sort
        samtools sort {params.basename}_nomerge.paired.unsorted.tmp.bam -o {output}

        # Then remove the intermediate, unsorted .bam
        rm {params.basename}_nomerge.paired.unsorted.tmp.bam
        """

## map_unmerged_unpaired: map trimmed, unmerged, unpaired reads to reference
#   BWA mem algorithm. Settings:
#   -M Mark shorter split hits as secondary (for Picard compatibility).
#   -t number of threads
#   -R read group, added through lambda function
#   then uses samtools view and samtools sort to convert to bam and sort
#   samtools view options:
#   -b output in bam format
rule map_unmerged_unpaired:
    input:
        reads_forward="processed_reads/trimmed/{sample}.nopair.R1.fq.qz",
        reads_reverse="processed_reads/trimmed/{sample}.nopair.R2.fq.qz",
        genome=REF
    output:
        mapped_forward = "processed_reads/mapped/{sample}.nopair.R1.sorted.bam",
        mapped_reverse = "processed_reads/mapped/{sample}.nopair.R2.sorted.bam"
    params:
        basename="processed_reads/mapped/{sample}",
        read_group=make_RG
    log:
        forward="logs/mapping/{sample}_nopair_R1.log",
        rev="logs/mapping/{sample}_nopair_R2.log"
    resources:
        cpus=8
    shell:
        """
        
        # Run bwa mem, pipe to samtools view to convert to bam, save as a tmp.bam
        # Read 1
        bwa mem -M -t {resources.cpus} -R '{params.read_group}' {input.genome} {input.reads_forward} 2> {log.forward} | samtools view -b - > {params.basename}_nopair.R1.unsorted.tmp.bam
        # Read 2
        bwa mem -M -t {resources.cpus} -R '{params.read_group}' {input.genome} {input.reads_reverse} 2> {log.rev} | samtools view -b - > {params.basename}_nopair.R2.unsorted.tmp.bam

        # then sort
        samtools sort {params.basename}_nopair.R1.unsorted.tmp.bam -o {output.mapped_forward}
        samtools sort {params.basename}_nopair.R2.unsorted.tmp.bam -o {output.mapped_reverse}

        # Then remove the intermediate, unsorted .bam
        rm {params.basename}_nopair.R1.unsorted.tmp.bam
        rm {params.basename}_nopair.R2.unsorted.tmp.bam
        """

## merge_bams_by_sample : merge bam files by sample and run
# merges bams across the 4 types of mapped reads (assembled, paired unassembled, and unpaired SEs)
# for a given sample/lane/sequencing run combination
# use samtools merge, -t is threads, rest is default
rule merge_sample_bams:
    input: 
        merged="processed_reads/mapped/{sample}.merged.sorted.bam",
        unmerged_pair="processed_reads/mapped/{sample}.nomerge.paired.sorted.bam",
        nopair_fwd="processed_reads/mapped/{sample}.nopair.R1.sorted.bam",
        nopair_rev="processed_reads/mapped/{sample}.nopair.R2.sorted.bam"
    log:
        "logs/merge_bams/{sample}_merge.log"
    resources:
        cpus=8
    output:
        "processed_reads/per_sample_bams/{sample}.sorted.bam"
    shell:
        """
        samtools merge -t {resources.cpus} {output} {input.merged} {input.unmerged_pair} {input.nopair_fwd} {input.nopair_rev} 2> {log}
        """

## index_deduped_bams: index bams
rule index_raw_bams:
    input:
        "processed_reads/per_sample_bams/{sample}.sorted.bam"
    output:
        "processed_reads/per_sample_bams/{sample}.sorted.bam.bai"
    resources:
        cpus=1
    shell:
        """
        samtools index -b {input}
        """


## remove_duplicates: remove duplicates with Picard
# Remove duplicates with Picard. options:
# REMOVE_DUPLICATES=true; remove duplicates, instead of default behavior (which outputs them and flags duplicated reads)
rule remove_duplicates:
    input:
        "processed_reads/per_sample_bams/{sample}.sorted.bam"
    output:
        bam="processed_reads/dedup/{sample}.sorted.bam",
        metrics="logs/remove_duplicates/{sample}_dedup_metrics.txt"
    log:
        log="logs/remove_duplicates/{sample}_dedup_log.log"
    resources:
        cpus=1
    shell:
        """
        picard MarkDuplicates INPUT={input} METRICS_FILE={output.metrics} OUTPUT={output.bam} REMOVE_DUPLICATES=true 2> {log.log}
        """


## multiqc_dedup_bam_report: collate qualimap reports on dedup bams
rule multiqc_dedup_report:
    input:
        expand("logs/remove_duplicates/{sample}_dedup_metrics.txt", sample = samples)
    output:
        "processed_reads/QC/multiqc/multiqc_report_dedup_Picard.html"
    shell:
        """
        multiqc logs/remove_duplicates/ -o processed_reads/QC/multiqc/ -n multiqc_report_dedup_Picard.html
        """
    



## index_deduped_bams: index bams
rule index_dedup_bams:
    input:
        "processed_reads/dedup/{sample}.sorted.bam"
    output:
        "processed_reads/dedup/{sample}.sorted.bam.bai"
    resources:
        cpus=1
    shell:
        """
        samtools index -b {input}
        """


## qualimap_final_bam: run qualimap on final bam file
# default options, only changed number of threads with -nt
rule qualimap_raw_bam:
    input:
        bam="processed_reads/per_sample_bams/{sample}.sorted.bam",
        bai="processed_reads/per_sample_bams/{sample}.sorted.bam.bai"
    output:
        "processed_reads/QC/qualimap/raw_bams/{sample}/qualimapReport.html"
    resources:
        cpus=8
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {resources.cpus} -outdir processed_reads/QC/qualimap/raw_bams/{wildcards.sample}/ -outformat html --java-mem-size=4G
        """

## qualimap_final_bam: run qualimap on final bam file
# default options, only changed number of threads with -nt
rule qualimap_dedup_bam:
    input:
        bam="processed_reads/dedup/{sample}.sorted.bam",
        bai="processed_reads/dedup/{sample}.sorted.bam.bai"
    output:
        "processed_reads/QC/qualimap/dedup_bams/{sample}/qualimapReport.html"
    resources:
        cpus=8
    shell:
        """
        qualimap bamqc -bam {input.bam} -nt {resources.cpus} -outdir processed_reads/QC/qualimap/dedup_bams/{wildcards.sample}/ -outformat html --java-mem-size=4G
        """

## multiqc_raw_bam_report: collate qualimap reports on raw bams
rule multiqc_raw_bam_report:
    input:
        expand("processed_reads/QC/qualimap/raw_bams/{sample}/qualimapReport.html", sample = samples)
    output:
        "processed_reads/QC/multiqc/multiqc_report_raw_bams.html"
    shell:
        """
        multiqc processed_reads/QC/qualimap/raw_bams -o processed_reads/QC/multiqc/ -n multiqc_report_raw_bams.html
        """
    
## multiqc_dedup_bam_report: collate qualimap reports on dedup bams
rule multiqc_dedup_bam_report:
    input:
        expand("processed_reads/QC/qualimap/dedup_bams/{sample}/qualimapReport.html", sample = samples)
    output:
        "processed_reads/QC/multiqc/multiqc_report_dedup_bams.html"
    shell:
        """
        multiqc processed_reads/QC/qualimap/dedup_bams -o processed_reads/QC/multiqc/ -n multiqc_report_dedup_bams.html
        """
    
