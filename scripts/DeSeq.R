#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(edgeR)
library(rtracklayer)
library(dplyr)


option_list <- list(
  make_option(c("-r", "--reads_table"), action = "store",
              help = "reads table generated from featureCounts"),
  make_option(c("-p", "--p_cutoff"), action = "store",
              help = "Pvalue cutoff", default = 0.05),
  make_option(c("-e", "--exp_desn"), action = "store",
              help = "experimental design file that contains information on groups"),
  make_option(c("-n", "--name"), action = "store",
              help = "name of feature from gff file that was chosen to represent each feature"),
  make_option(c("-o", "--out_dir"), action = "store",
              help = "an output directory where all outputs will be stored")
)

opt <- parse_args(OptionParser(option_list = option_list))

reads_file <- opt$reads_table
pcutoff <- opt$p_cutoff
group_file <- opt$exp_desn
out_dir <- opt$out_dir

# read the output of featureCounts
read.counts <- read.table(reads_file, sep = "\t", header = TRUE)

# rename the column headers
names(read.counts) <- gsub(".*mapping_results.", "", names(read.counts),
                           perl = TRUE)
names(read.counts) <- gsub("_srt.bam", "", names(read.counts), perl = TRUE)

# assign row names as gene names
row.names(read.counts) <- read.counts[, 1]

# gene information
gene.info <- read.counts[, c(1:6)]

# remove first six columns that have metadata
read.counts <- read.counts[, -c(1:6)]


# read in the table with group info
group_table <- read.delim(group_file, row.names = 1)
group_table <- select(group_table, Group)

read.counts <- read.counts[, rownames(group_table)]
deseq_ds <- DESeq2::DESeqDataSetFromMatrix(countData = read.counts,
                                           colData = group_table,
                                           design = ~ Group)

# remove genes without any counts
deseq_ds <- deseq_ds[rowSums(counts(deseq_ds)) > 0, ]

# calculate size factors
dds <- DESeq2::DESeq(deseq_ds)

# create the output directory
ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))

# PCA plot of samples
pdf(file.path(out_dir, "PCA.pdf"))
DESeq2::plotPCA(DESeq2::varianceStabilizingTransformation(dds),
                ntop = 1000, intgroup = c("Group"))
dev.off()



# create heatmap of count data.
pdf(file.path(out_dir, "count_heatmap.pdf"))
heatmap(as.matrix(read.counts))
dev.off()

# convert to a DGElist, remove rows where there are no reads mapped
edger_dge <- edgeR::DGEList(counts = read.counts, group = group_table$Group,
                            remove.zeros = TRUE, genes = gene.info)

# histogram of count per million
pdf(file.path(out_dir, "cpm_histograpm.pdf"))
hist(log2(rowSums(cpm(edger_dge))))
dev.off()

# only keep genes/features that have total cpm in at least all samples
keep <- rowSums(cpm(edger_dge) > 1) >= nrow(group_table)
edger_dge <- edger_dge[keep, ]

# recompute library size after filtering
edger_dge$samples$lib.size <- colSums(edger_dge$counts)

# calculate the normalization factors for the library size
edger_dge <- edgeR::calcNormFactors(edger_dge, method = "TMM")

# estimate the dispersion for all read counts across all samples
edger_dge <- edgeR::estimateDisp(edger_dge, design.robust = TRUE)

# calculate RPKM and CPM
rpkm_results <- edgeR::rpkm(edger_dge)
cpm_results <- edgeR::cpm(edger_dge)

out_rpkm <- file.path(out_dir, "RPKM.csv")
out_cpm <- file.path(out_dir, "CPM.csv")

write.csv(rpkm_results, file = out_rpkm)
write.csv(cpm_results, file = out_cpm)


# Create pairwise comparisons to find DGEs
pair.comb <- function(exp_des){
    # get all pariwise combination from experimental design file
    exp_desn <- read.table(exp_des, sep = "\t", header = TRUE)
    categories <- unique(exp_desn$Group)
    pairs <- combn(categories, 2, simplify = FALSE)
    return(pairs)
}

# get all possible pairs
all_pairs <- pair.comb(group_file)

for (n in 1:length(all_pairs) ) {
  # filename strings for each comparisons
  filename <- paste(all_pairs[[n]][1], all_pairs[[n]][2], "et.csv", sep = "__")
  filename_sig <- paste(all_pairs[[n]][1], all_pairs[[n]][2], "sig.csv", sep = "__")

  # sample matrix
  pair1 <- as.character(all_pairs[[n]][1])
  pair2 <- as.character(all_pairs[[n]][2])
  pairs <- c(pair1, pair2)

  # exact test
  deseq_diff <- DESeq2::results(dds, contrast = c("Group", pair1, pair2))
  deseq_diff <- deseq_diff[order(as.numeric(deseq_diff$pvalue)),]
  deseq_sig <- subset(deseq_diff, pvalue < as.numeric(pcutoff))

  # full path to output directory
  out_table <- file.path(out_dir, filename)
  out_table_sig <- file.path(out_dir, filename_sig)

  # write the file
  write.csv(deseq_diff, out_table)
  write.csv(deseq_sig, out_table_sig)

}