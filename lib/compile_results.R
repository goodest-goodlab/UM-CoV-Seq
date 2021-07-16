# Load packages -----------------------------------------------------------
library(dplyr)
library(stringr)
library(tidyr)

# Get arguments -----------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# 1 consensus stat
# 2 multiqc results
# 3 variant directory
# 4 pangolin results
# 5 nextclade results
# 6 output file
stats_file <- args[1]
multiqc_file <- args[2]
variant_dir <- args[3]
pangolin_file <- args[4]
nextclade_file <- args[5]
vcf_dir <- args[6]
output_file <- args[7]

# Process the stats file
stats <- read.csv(stats_file)  %>% 
  separate(sample, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, N, perc_N, n_less5k)

# Process the coverage file
coverage <- read.delim(multiqc_file) %>% 
  select(sample = Sample, perc_ref_1X_cov = QualiMap_mqc.generalstats.qualimap.1_x_pc, perc_ref_30X_cov = QualiMap_mqc.generalstats.qualimap.30_x_pc)

# Count variants
tsvs <- list.files(variant_dir, pattern = ".tsv", full.names = T)
variants <- data.frame(sample = rep(NA, length(tsvs)), ivar_vars = rep(NA, length(tsvs)))
i <- 1
for (file in tsvs) {
  sample <- str_remove(basename(file), ".tsv")
  var_number <- read.delim(file, sep = "\t") %>% 
    group_by(POS) %>% 
    tally() %>% 
    ungroup() %>% 
    tally() %>% 
    pull(n)
  
  variants$sample[i] <- sample
  variants$ivar_vars[i] <- var_number
  i <- i + 1
}


# Pangolin
pangolin <- read.csv(pangolin_file) %>% 
  separate(taxon, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, pango_qc = status, pango_lineage = lineage) 

# Nextclade
nextclade <- read.delim(nextclade_file, sep = "\t", na.strings = "") %>% 
  separate(seqName, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, next_qc = qc.overallStatus, next_clade = clade) 


# Combine and output
stats %>% 
  left_join(coverage) %>%
  left_join(variants) %>% 
  left_join(pangolin) %>% 
  left_join(nextclade) %>% 
  write.csv(., file = output_file, row.names = F)
