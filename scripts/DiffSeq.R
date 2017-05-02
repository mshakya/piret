#!/usr/bin/env Rscript
#
# usage: Rscript ../../bin/DiffSeq.R -r sense_reads.csv -g IV_A_PR8_34_mofied_cds.gff -s sum_exp_stats.txt -e ../2016-10-14/ExpDescrNonPolyA_nogz.txt -t region -n segment
#


library(optparse)
library(edgeR)
library(rtracklayer)
library(dplyr)


option_list <- list(
  make_option(c("-r", "--reads_table"), action = "store",
              help = "reads table generated from featureCounts"),
  make_option(c("-g", "--gtf_file"), action = "store",
              help = "gtf file to extract information on length"),
  make_option(c("-e", "--exp_desn"), action = "store",
              help = "experimental design file that contains information on groups"),
  make_option(c("-t", "--type"), action = "store",
              help = "type of (gene/CDS/region) that was used to map"),
  make_option(c("-n", "--name"), action = "store",
              help = "name of feature from gff file that was chosen to represent each feature"),
  make_option(c("-d", "--tot_reads"), action = "store", default = "Proper_paired_prokaryote_reads",
              help = "header from sum_exp_stats.txt to be used for total reads mapped for RPKM"),
  make_option(c("-o", "--out_table"), action = "store",
              help = "output table (csv) with RPKM")
)

opt <- parse_args(OptionParser(option_list=option_list))

reads_file <- opt$reads_table
gtf <- opt$gtf_file
summ <- opt$summ_table
group_file <- opt$exp_desn
feat_name <- opt$name
tot_reads <- opt$tot_reads
feat_type <- opt$type
out_table <- opt$out_table


get_samp_name <- function(col_name) {
  # Function to extract sample name from the full path.
  if (grepl("/", col_name)) {
    samp_path <- strsplit(col_name, "/")[[1]]
    new_col <- samp_path[length(samp_path) - 2]
  return(new_col)
  } else
  {
    return(col_name)
    }
}

# read experimental design file
group_table <- read.delim(group_file)
group_table <- select(group_table, ID, Group)

# read count table
reads_table <- read.table(reads_file,
  check.names = FALSE, stringsAsFactors = FALSE, sep = "\t",
  header = TRUE)

# change the column name to just have sample names
col_list <- colnames(reads_table)
colnames(reads_table)<-sapply(col_list, get_samp_name)


# convert to DGE object
DGE_object <- DGEList(counts = reads_table[, 7:12], genes = reads_table[, 1:6],
                      group = group_table$Group)

pdf("~/Desktop/test.pdf")
plotMDS( DGE_object,  main = "MDS Plot for Count Data", labels = colnames( DGE_object$counts ) )
plotSmear(DGE_object)
dev.off()

