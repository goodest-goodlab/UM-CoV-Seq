# Load packages -----------------------------------------------------------

## First specify the packages of interest
packages = c("ggplot2", "ggrepel", "dplyr", "stringr", "tidyr")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)

# Functions -----------------------------------------------------------

bartheme <- function () {  
  # My standard theme for most plots.
  theme_classic() %+replace% 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=16), 
          axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black",angle=90), 
          axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0),color="black"),
          axis.line=element_line(colour='#595959',size=0.75),
          axis.ticks=element_line(colour="#595959",size = 1),
          axis.ticks.length=unit(0.2,"cm"),
          legend.position="right",
          legend.key.width = unit(0.75,  unit = "cm"),
          legend.spacing.x = unit(0.25, 'cm'),
          legend.title = element_blank(),
          legend.text=element_text(size=12),
          plot.title = element_text(hjust=0.5, size=20),
          plot.margin = unit(c(1,1,1,1), "cm")
    )
}

# Get arguments -----------------------------------------------------------

cat("\n")
cat(as.character(Sys.time()), " | Reading arguments\n")
args = commandArgs(trailingOnly=TRUE)

# Previously, loaded in all the files and
# base directories as arguments and parameters.
# Makes for a lot of argument calling and parsing. But, perhaps
# easier maintenance: if our file structure changes,
# can change it in the snakefile without editing the R script 
# (but we still need to edit the snakemake rule and arguments and make sure they're in order)
# I think I'll switch to a hopefully easier option: just pass the output file
# to R, and from there we can learn all that we need to know about the
# directory structure.
# only other argument is the masking flag, which tells us whether masking occured or not

output_file = args[1]
mask_problematic = args[2]

if (mask_problematic %in% c("TRUE", "FALSE") == F) {
  stop(paste0("mask_problematic flag is ", mask_problematic, ", not TRUE or FALSE"))
}

# From the ouput file, figure out all the other filenames we need:

# The base output directory
base_dir = dirname(output_file)
# The batch name to use
batch_name = str_remove(basename(output_file), "-summary.csv")

# The multiqc file with coverage info
multiqc_file = file.path(base_dir, 
                         "results",
                         "multiqc", 
                         "multiqc_report_raw_bams_data", 
                         "multiqc_general_stats.txt")

# Directories with variant calls for ivar and GATK
ivar_dir = file.path(base_dir, 
                     "results", 
                     "ivar")
# With GATK, the .vcf file we use depends on whether 
# problematic site masking was performed
if (mask_problematic == "TRUE") {
  gatk_dir = file.path(base_dir, 
                       "results",
                       "gatk", 
                       "vcfs",
                       "filtered_and_masked")
} else {
  gatk_dir = file.path(base_dir, 
                       "results",
                       "gatk",  
                       "vcfs",
                       "filtered")
}

# Stats files of consensus seqs
ivar_stats_file = file.path(base_dir, 
                       "results", 
                       "ivar",
                       "all_samples_consensus_stats.csv")
gatk_stats_file = file.path(base_dir, 
                       "results", 
                       "gatk",
                       "all_samples_consensus_stats.csv")


# Pangolin lineage reports
ivar_pangolin_file = file.path(base_dir, 
                          "results", 
                          "ivar-pangolin",
                          "lineage_report.csv")
gatk_pangolin_file = file.path(base_dir, 
                          "results", 
                          "gatk-pangolin",
                          "lineage_report.csv")


# The nextclade lineage reports
ivar_nextclade_file = file.path(base_dir, 
                          "results", 
                          "ivar-nextclade",
                          "nextclade_report.tsv")
gatk_nextclade_file = file.path(base_dir, 
                          "results", 
                          "gatk-nextclade",
                          "nextclade_report.tsv")

cat(as.character(Sys.time()), " | Output file:", output_file, "\n")
cat(as.character(Sys.time()), " | Reading sample IDs from file:", multiqc_file, "\n")


samples = read.delim(multiqc_file) %>% 
  select(sample = Sample)

###############

cat(as.character(Sys.time()), " | Getting coverage stats\n")
# Process the coverage file
coverage <- read.delim(multiqc_file) %>% 
  select(sample = Sample, perc_ref_1X_cov = QualiMap_mqc.generalstats.qualimap.1_x_pc, 
         perc_ref_30X_cov = QualiMap_mqc.generalstats.qualimap.30_x_pc, 
         avg.coverage = QualiMap_mqc.generalstats.qualimap.mean_coverage,
         perc.reads.mapped = QualiMap_mqc.generalstats.qualimap.percentage_aligned)
# Read the coverage file and select columns of interest

coverage$cov.filter = ifelse(coverage$perc_ref_1X_cov < 90, "FILTER", "PASS")
# Add a column for the coverage at 1X filter

cat(as.character(Sys.time()), " | Generating 1x coverage plot\n")
cov_perc_1x_p = ggplot(coverage, aes(x=sample, y=perc_ref_1X_cov)) +
  geom_segment(aes(x=sample, y=100, xend=sample, yend=perc_ref_1X_cov), linetype="dotted", color="#666666") +
  geom_point(size=2, color="#920000") +
  geom_hline(yintercept=90, size=1, linetype="dashed", color="#56b4e9") +
  geom_text_repel(aes(label=ifelse(perc_ref_1X_cov < 90, as.character(sample), "")), show.legend=FALSE) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_discrete(position = "top") +
  ylab("% of sites covered\nby at least 1 read") +
  xlab("sample") +
  ggtitle(paste0("Analysis batch: ", batch_name)) + 
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=-0.001, size=6))
# Plot the percent of sites that have at least 1 read covering them

cov_perc_1x_p_out = file.path(base_dir, 
                              paste0(batch_name,"-1x-cov.pdf"))  
cat(as.character(Sys.time()), " | Saving 1x coverage plot:", cov_perc_1x_p_out, "\n")
ggsave(cov_perc_1x_p_out, cov_perc_1x_p, width=10, height=4, units="in")
# Save the plot

cat(as.character(Sys.time()), " | Generating average coverage plot\n")
coverage$sample = factor(coverage$sample, levels=coverage$sample[order(coverage$avg.coverage, decreasing=T)])
# This orders the columns by read depth in descending order
avg_cov_p = ggplot(coverage, aes(x=sample, y=avg.coverage)) +
  geom_bar(stat="identity", fill="transparent", color="#333333") +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept=mean(coverage$avg.coverage, na.rm=T), linetype="dashed", color="#999999") +
  geom_hline(yintercept=median(coverage$avg.coverage, na.rm=T), linetype="dashed", color="#920000") +
  xlab("sample") +
  ggtitle(paste0("Analysis batch: ", batch_name)) + 
  ylab("Avg. depth") +
  bartheme() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=6))
# Plot the average read depth per sample


avg_cov_out = file.path(base_dir, 
                        paste0(batch_name, "-avg-cov.pdf"))  
cat(as.character(Sys.time()), " | Saving average coverage plot:", avg_cov_out, "\n")
ggsave(avg_cov_out, avg_cov_p, width=10, height=4, units="in")
# Save the plot

##

cat(as.character(Sys.time()), " | Generating % reads mapped plot\n")
perc_reads_p = ggplot(coverage, aes(x=sample, y=perc.reads.mapped)) +
  geom_segment(aes(x=sample, y=100, xend=sample, yend=perc.reads.mapped), linetype="dotted", color="#666666") +
  geom_point(size=2, color="#920000") +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_discrete(position = "top") +
  ylab("% reads mapped") +
  xlab("sample") +
  ggtitle(paste0("Analysis batch: ", batch_name)) + 
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=-0.001, size=6))
#print(perc_reads_p)
# Plot the percent of reads mapped per sample

perc_reads_p_out = file.path(base_dir, 
                             paste0(batch_name, "-percent-mapped.pdf"))  
cat(as.character(Sys.time()), " | Saving % reads mapped plot:", perc_reads_p_out, "\n")
ggsave(perc_reads_p_out, perc_reads_p, width=10, height=4, units="in")
# Save the plot

###############
# Count variants

# Count ivar variants
cat(as.character(Sys.time()), " | Counting ivar variants\n")
tsvs <- list.files(ivar_dir, pattern = ".tsv", full.names = T)
ivar_vars = data.frame("sample"=c(), "ivar.num.vars"=c(), "ivar.vars"=c())
i <- 1
for (file in tsvs) {
  sample <- str_remove(basename(file), ".tsv")
  cur_vars_df = read.delim(file, sep = "\t") %>% group_by(POS) %>% slice(1)
  cur_vars = c()
  for(j in 1:nrow(cur_vars_df)){
    cur_var = paste(cur_vars_df[j,]$POS, cur_vars_df[j,]$ALT, sep=":")
    cur_vars = c(cur_vars, cur_var)
  }
  ivar_vars = rbind(ivar_vars, data.frame("sample"=sample, 
                                          "ivar.num.vars"=nrow(cur_vars_df), 
                                          "ivar.vars"=I(list(cur_vars))))
  i <- i + 1
}


cat(as.character(Sys.time()), " | Counting gatk variants\n")
vcfs = list.files(gatk_dir, pattern=".fvcf.gz$", full.names=T)
gatk_vars = data.frame("sample"=c(), "gatk.num.vars"=c(), "gatk.vars"=c())
for(file in vcfs){
  sample = str_remove(str_remove(basename(file), ".fvcf.gz"), ".masked")
  cur_vars_df = read.csv(file, comment.char="#", sep="\t", header=F)
  cur_vars_df = subset(cur_vars_df, V5!= "." & V5 != "*" & V5 != "N" & V7 == "PASS")
  cur_vars = c()
  for(j in 1:nrow(cur_vars_df)){
    cur_var = paste(cur_vars_df[j,]$V2, cur_vars_df[j,]$V5, sep=":")
    cur_vars = c(cur_vars, cur_var)
  }
  gatk_vars = rbind(gatk_vars, data.frame("sample"=sample, 
                                          "gatk.num.vars"=nrow(cur_vars_df), 
                                          "gatk.vars"=I(list(cur_vars))))
  i <- i + 1
}
# Count gatk variants

cat(as.character(Sys.time()), " | Combining and comparing variant calls\n")
variants = merge(ivar_vars, gatk_vars, by="sample")
variants$num.shared = NA
# Merge ivar and gatk variant counts

for(i in 1:nrow(variants)){
  cur_sample = variants[i,]
  #print(cur_sample$sample)
  ivar_vars = data.frame("sample"=cur_sample$sample, "var"=unlist(cur_sample$ivar.vars))
  gatk_vars = data.frame("sample"=cur_sample$sample, "var"=unlist(cur_sample$gatk.vars))
  
  ivar_uniq = ivar_vars[!ivar_vars$var %in% gatk_vars$var,]
  gatk_uniq = gatk_vars[!gatk_vars$var %in% ivar_vars$var,]
  shared = ivar_vars[ivar_vars$var %in% gatk_vars$var,]

  variants[i,]$num.shared = nrow(shared)
}
# Count the number of variants that are shared between ivar and gatk for all samples
# May want to plot this somehow? But with so much agreement that seems pointless right now

###############
# Process the consensus stats files

cat(as.character(Sys.time()), " | Getting ivar consensus stats\n")
ivar_stats <- read.csv(ivar_stats_file, header=T) %>% 
  select(sample, N, length, perc.n, n.filter)
names(ivar_stats) = c("sample", "ivar.Ns", "ivar.length", "ivar.perc.n", "ivar.n.filter")
# Read the ivar stats

cat(as.character(Sys.time()), " | Getting gatk consensus stats\n")
gatk_stats <- read.csv(gatk_stats_file, header=T) %>% 
  select(sample, N, length, perc.n, n.filter)
names(gatk_stats) = c("sample", "gatk.Ns", "gatk.length", "gatk.perc.n", "gatk.n.filter")
# Read the gatk stats

cat(as.character(Sys.time()), " | Combining consensus stats\n")
stats = merge(ivar_stats, gatk_stats, by="sample")
# Merge the ivar and gatk consensus stats

cat(as.character(Sys.time()), " | Generating Ns plot\n")
n_p = ggplot(stats, aes(x=sample, y=ivar.Ns, color="iVar")) +
  geom_segment(aes(x=sample, y=0, xend=sample, yend=ivar.Ns), linetype="dotted", color="#666666") +
  geom_point(size=2, color="#920000") +
  geom_segment(aes(x=sample, y=0, xend=sample, yend=gatk.Ns), linetype="dotted", color="#666666") +
  geom_point(aes(x=sample, y=gatk.Ns, color="GATK"), size=2) +
  geom_hline(yintercept=5000, size=1, linetype="dashed", color="#56b4e9") +
  geom_text_repel(aes(label=ifelse(ivar.n.filter=="FILTER", as.character(sample), "")), show.legend=FALSE) +
  geom_text_repel(aes(label=ifelse(gatk.n.filter=="FILTER", as.character(sample), "")), show.legend=FALSE) +
  scale_y_continuous(expand=c(0,0), limits=c(0,31000)) +
  ylab("# Ns") +
  xlab("sample") +
  ggtitle(paste0("Analysis batch: ", batch_name)) + 
  scale_color_manual(name="", values=c("iVar"="#920000", "GATK"="#db6d00")) +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=6),
        legend.position="bottom")
#print(perc_n_p)
# Plot the number of Ns for each sample

n_p_out = file.path(base_dir, 
                    paste0(batch_name, "-num-Ns.pdf")) 
cat(as.character(Sys.time()), " | Saving Ns plot:", n_p_out, "\n")
ggsave(n_p_out, n_p, width=10, height=4, units="in")
# Save the plot

###############
# Pangolin

cat(as.character(Sys.time()), " | Reading ivar pangolin lineages\n")
ivar_pangolin <- read.csv(ivar_pangolin_file) %>% 
  separate(taxon, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, ivar_pango_qc = status, ivar_pango_lineage = lineage) 
# Read the ivar lineages

cat(as.character(Sys.time()), " | Reading gatk pangolin lineages\n")
gatk_pangolin <- read.csv(gatk_pangolin_file) %>% 
  separate(taxon, into = c("sample", "ext1", "ext2", "ext3", "ext4", "ext5", "ext6",
                           "ext7", "ext8", "ext9", "ext10", "ext11", "ext12", "ext13"), sep = "_") %>% 
  select(sample, gatk_pango_qc = status, gatk_pango_lineage = lineage) 
# Read the gatk lineages

cat(as.character(Sys.time()), " | Combining and comparing pangolin lineages\n")
pangolin = merge(ivar_pangolin, gatk_pangolin, by="sample")
pangolin$pango.match = ifelse(as.character(pangolin$gatk_pango_lineage) == as.character(pangolin$ivar_pango_lineage), 1, 0)
# Merge the ivar and gatk tables and add a column based on whether they match or not

###############
# Nextclade

cat(as.character(Sys.time()), " | Reading ivar nextclade lineages\n")
ivar_nextclade <- read.delim(ivar_nextclade_file, sep = "\t", na.strings = "") %>% 
  separate(seqName, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, ivar_next_qc = qc.overallStatus, ivar_next_clade = clade) 
# Read the ivar lineages

cat(as.character(Sys.time()), " | Reading gatk nextclade lineages\n")
gatk_nextclade <- read.delim(gatk_nextclade_file, sep = "\t", na.strings = "") %>% 
  separate(seqName, into = c("sample", "ext1", "ext2"), sep = "_") %>% 
  select(sample, gatk_next_qc = qc.overallStatus, gatk_next_clade = clade) 
# Read the gatk lineages

cat(as.character(Sys.time()), " | Combining and comparing nextclade lineages\n")
nextclade = merge(ivar_nextclade, gatk_nextclade, by="sample")
nextclade$nextclade.match = ifelse(as.character(nextclade$gatk_next_clade) == as.character(nextclade$ivar_next_clade), 1, 0)
# Merge the ivar and gatk tables and add a column based on whether they match or not

###############
# Combine tables
# Tables at this point should includes: stats, coverage, variants, pangolin, nextclade

variants = subset(variants, select=-c(ivar.vars, gatk.vars))
# Remove the list of variants from the variants table before merging to final table

cat(as.character(Sys.time()), " | Combining all tables\n")
final_table = merge(coverage, stats, by="sample")
final_table = merge(final_table, variants, by="sample")
final_table = merge(final_table, pangolin, by="sample")
final_table = merge(final_table, nextclade, by="sample")

# Merge all the tables

cat(as.character(Sys.time()), " | Writing summary output file:", output_file, "\n")
write.csv(final_table, file=output_file, row.names=F)
# Write the summary output file

cat("------------------------------------\n")
passing_samples = subset(final_table, ivar.n.filter=="PASS" & gatk.n.filter=="PASS" & cov.filter == "PASS")
cat(as.character(Sys.time()), " | NUMBER OF SAMPLES THAT PASS ALL FILTERS IN BOTH PIPELINES:", nrow(passing_samples), "/", nrow(final_table), "\n")
cat(as.character(Sys.time()), " | PASSING SAMPLES WITH PANGOLIN LINEAGE MATCHES:            ", sum(passing_samples$pango.match), "/", nrow(passing_samples), "\n")
cat(as.character(Sys.time()), " | PASSING SAMPLES WITH NEXTCLADE LINEAGE MATCHES:           ", sum(passing_samples$nextclade.match), "/", nrow(passing_samples), "\n")
# Print out some counts about samples that pass filters and matching lineages.

###############

