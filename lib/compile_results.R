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

sample_file = args[1]
batch = args[2]
# Overall sample file and batch id

batch_dir = args[3]
# The base batch output directory (should be results/BATCH)

multiqc_file = args[4]
# The multiqc file with coverage info

ivar_dir = args[5]
# The directories containing the variant calls 

ivar_cons_file = args[6]
# The stats files for the consensus sequences

ivar_pangolin_file = args[7]

###############
# Determin output file

output_file = paste(batch_dir, batch, "-summary.csv", sep="")
cat(as.character(Sys.time()), " | Output file:", output_file, "\n")

###############
# Read the sample file

cat(as.character(Sys.time()), " | Reading samples file:", sample_file, "\n")
samples = read.csv(sample_file, comment.char="#", quote='"')
samples = subset(samples, MiSeq.Run == batch)

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
  #geom_text_repel(aes(label=sample)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_discrete(position = "top") +
  ylab("% of sites covered\nby at least 1 read") +
  xlab(paste(batch, "sample")) +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=-0.001, size=6))
#print(cov_perc_1x_p)
# Plot the percent of sites that have at least 1 read covering them

cov_perc_1x_p_out = paste(batch_dir, batch, "-1x-cov.pdf", sep="")
cat(as.character(Sys.time()), " | Saving 1x coverage plot:", cov_perc_1x_p_out, "\n")
ggsave(cov_perc_1x_p_out, cov_perc_1x_p, width=10, height=4, units="in")
# Save the plot

##

cat(as.character(Sys.time()), " | Generating average coverage plot\n")
coverage$sample = factor(coverage$sample, levels=coverage$sample[order(coverage$avg.coverage, decreasing=T)])
# This orders the columns by read depth in descending order
avg_cov_p = ggplot(coverage, aes(x=sample, y=avg.coverage)) +
  geom_bar(stat="identity", fill="transparent", color="#333333") +
  #scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_hline(yintercept=mean(coverage$avg.coverage, na.rm=T), linetype="dashed", color="#999999") +
  geom_hline(yintercept=median(coverage$avg.coverage, na.rm=T), linetype="dashed", color="#920000") +
  xlab(paste(batch, "sample")) +
  ylab("Avg. depth") +
  bartheme() + 
  theme(axis.text.x=element_text(angle=45, hjust=1, size=6))
#print(avg_cov_p)
# Plot the average read depth per sample

avg_cov_out = paste(batch_dir, batch, "-avg-cov.pdf", sep="")
cat(as.character(Sys.time()), " | Saving average coverage plot:", avg_cov_out, "\n")
ggsave(avg_cov_out, avg_cov_p, width=10, height=4, units="in")
# Save the plot

##

cat(as.character(Sys.time()), " | Generating % reads mapped plot\n")
perc_reads_p = ggplot(coverage, aes(x=sample, y=perc.reads.mapped)) +
  geom_segment(aes(x=sample, y=100, xend=sample, yend=perc.reads.mapped), linetype="dotted", color="#666666") +
  geom_point(size=2, color="#920000") +
  #geom_hline(yintercept=5000, size=1, linetype="dashed", color="#56b4e9") +
  scale_y_continuous(expand=c(0,0), limits=c(0,100)) +
  scale_x_discrete(position = "top") +
  ylab("% reads mapped") +
  xlab(paste(batch, "sample")) +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=-0.001, size=6))
#print(perc_reads_p)
# Plot the percent of reads mapped per sample

perc_reads_p_out = paste(batch_dir, batch, "-percent-mapped.pdf", sep="")
cat(as.character(Sys.time()), " | Saving % reads mapped plot:", perc_reads_p_out, "\n")
ggsave(perc_reads_p_out, perc_reads_p, width=10, height=4, units="in")
# Save the plot

###############
# Count variants

cat(as.character(Sys.time()), " | Counting ivar variants\n")
tsvs <- list.files(ivar_dir, pattern = ".tsv", full.names = T)
#variants <- data.frame(sample = rep(NA, length(tsvs)), ivar_vars = rep(NA, length(tsvs)))
variants = data.frame("sample"=c(), "ivar.num.vars"=c())#, "ivar.vars"=c())
i <- 1
for (file in tsvs) {
  sample <- str_remove(basename(file), ".tsv")
  cur_vars_df = read.delim(file, sep = "\t") %>% group_by(POS) %>% slice(1)
  cur_vars = c()
  for(j in 1:nrow(cur_vars_df)){
    cur_var = paste(cur_vars_df[j,]$POS, cur_vars_df[j,]$ALT, sep=":")
    cur_vars = c(cur_vars, cur_var)
  }
  variants = rbind(variants, data.frame("sample"=sample, 
                                          "ivar.num.vars"=nrow(cur_vars_df)))
  i <- i + 1
}
# Count ivar variants

cat(as.character(Sys.time()), " | Getting ivar consensus stats\n")
stats <- read.csv(ivar_cons_file, header=T) %>% 
  select(sample, N, length, perc.n, n.filter)
names(stats) = c("sample", "ivar.Ns", "ivar.length", "ivar.perc.n", "ivar.n.filter")

cat(as.character(Sys.time()), " | Generating Ns plot\n")
n_p = ggplot(stats, aes(x=sample, y=ivar.Ns, color="iVar")) +
  geom_segment(aes(x=sample, y=0, xend=sample, yend=ivar.Ns), linetype="dotted", color="#666666") +
  geom_point(size=2, color="#920000") +
  geom_hline(yintercept=5000, size=1, linetype="dashed", color="#56b4e9") +
  geom_text_repel(aes(label=ifelse(ivar.n.filter=="FILTER", as.character(sample), "")), show.legend=FALSE) +
  scale_y_continuous(expand=c(0,0), limits=c(0,31000)) +
  ylab("# Ns") +
  xlab(paste(batch, "sample")) +
  scale_color_manual(name="", values=c("iVar"="#920000")) +
  bartheme() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=6),
        legend.position="bottom")
#print(perc_n_p)
# Plot the number of Ns for each sample

n_p_out = paste(batch_dir, batch, "-num-Ns.pdf", sep="")
cat(as.character(Sys.time()), " | Saving Ns plot:", n_p_out, "\n")
ggsave(n_p_out, n_p, width=10, height=4, units="in")
# Save the plot

###############
# Pangolin

cat(as.character(Sys.time()), " | Reading ivar pangolin lineages\n")
pangolin <- read.csv(ivar_pangolin_file) %>% 
  separate(taxon, into = c("ext", "sample", "ext2", "ext3", "ext4", "ext5"), sep = "_") %>% 
  select(sample, ivar_pango_qc = status, ivar_pango_lineage = lineage) 
# Read the ivar lineages


###############
# Combine tables
# Tables at this point should includes: stats, coverage, variants, pangolin, nextclade

#variants = subset(variants, select=-c(ivar.vars, gatk.vars))
# Remove the list of variants from the variants table before merging to final table

cat(as.character(Sys.time()), " | Combining all tables\n")
final_table = merge(coverage, stats, by="sample")
final_table = merge(final_table, variants, by="sample")
final_table = merge(final_table, pangolin, by="sample")
final_table$sample = str_replace(final_table$sample, "-ICV", "")
# Merge all the tables

cat(as.character(Sys.time()), " | Writing summary output file:", output_file, "\n")
write.csv(final_table, file=output_file, row.names=F)
# Write the summary output file

cat("------------------------------------\n")
passing_samples = subset(final_table, ivar.n.filter=="PASS" & cov.filter == "PASS")
cat(as.character(Sys.time()), " | NUMBER OF SAMPLES THAT PASS ALL FILTERS IN BOTH PIPELINES:", nrow(passing_samples), "/", nrow(final_table), "\n")
# Print out some counts about samples that pass filters and matching lineages.

###############

