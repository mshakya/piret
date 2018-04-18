#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
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
feature_name <- opt$name
out_dir <- opt$out_dir

# create the output directory
ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))

# read the output of featureCounts
read.counts <- read.table(reads_file, sep = "\t", header=TRUE, row.names=1)

# # rename the column headers
names(read.counts) <- gsub(".*mapping_results.", "", names(read.counts),
                           perl = TRUE)
names(read.counts) <- gsub("_srt.bam", "", names(read.counts), perl = TRUE)

# # assign row names as gene names
# row.names(read.counts) <- read.counts[, 1]


# gene information
gene.info <- read.counts[, c(1:5)]
print(head(gene.info))
gene.ranges <- makeGRangesFromDataFrame(gene.info)
print(gene.ranges)
# # remove first six columns that have metadata
# read.counts <- read.counts[, -c(1:6)]

# # read in the table with group info
group_table <- read.delim(group_file, row.names = 1)
group_table <- select(group_table, Group)

read.counts <- read.counts[, rownames(group_table)]
deseq_ds <- DESeq2::DESeqDataSetFromMatrix(countData = read.counts,
                                           colData = group_table,
                                           design = ~ Group,
                                           tidy = FALSE,
                                           rowRanges=gene.ranges)

# remove genes without any counts
deseq_ds <- deseq_ds[rowSums(counts(deseq_ds)) > 0, ]

# calculate size factors
dds <- DESeq2::DESeq(deseq_ds)

#TODO: Need to make sure that all gff lines have IDs
# ############### calculate FPKM and FPM #########################################
fpkm_results <- DESeq2::fpkm(deseq_ds, robust=TRUE)
fpm_results <- DESeq2::fpm(deseq_ds)
out_fpkm <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_FPKM.csv", sep=""))
out_fpm <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_FPM.csv", sep=""))
write.csv(fpkm_results, file = out_fpkm)
write.csv(fpm_results, file = out_fpm)
################################################################################


################ histogram of count per million ################################
out_fpm_hist <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_fpm_histogram.pdf", sep=""))
fpm_results <- dplyr::filter_all(as.data.frame(fpm_results), any_vars(. != 0))
fpm_data <- reshape2::melt(as.data.frame(fpm_results), variable.name="sample", value.name="CPM")
fpm_hist <- ggplot(data = fpm_data, mapping = aes(x = CPM)) +  theme_bw() +
        geom_histogram(bins=100) + xlab("CPM") + ylab(feature_name) + facet_wrap(~sample)
ggsave(out_fpm_hist, fpm_hist, device = "pdf")

################## boxplot of fragments per million ################################
group_table2 <- add_rownames(group_table, "sample")
fpm_data_boxplot <- merge(x=fpm_data, y=group_table2)
out_fpm_violin <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_fpm_violin.pdf", sep=""))
fpm_violin <- ggplot(data = fpm_data_boxplot, mapping = aes(x=sample, y=CPM)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
fpm_violin_group <- ggplot(data = fpm_data_boxplot, mapping = aes(x=Group, y=CPM)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
pdf(out_fpm_violin)
fpm_violin
fpm_violin_group
dev.off()

##############################histogram of fpkm ################################
out_fpkm_hist <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_fpkm_histogram.pdf", sep=""))
fpkm_results <- dplyr::filter_all(as.data.frame(fpkm_results), any_vars(. != 0))
fpkm_data <- reshape2::melt(as.data.frame(fpkm_results), variable.name="sample", value.name="fpkm")
fpkm_hist <- ggplot(data=fpkm_data, mapping=aes(x=fpkm)) +  theme_bw() +
        geom_histogram(bins=100) + xlab("fpkm") + ylab(feature_name) + facet_wrap(~sample)
ggsave(out_fpkm_hist, fpkm_hist, device = "pdf")

##############################heatmap of fpkm ################################
out_fpkm_heatmap <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_fpkm_heatmap.pdf", sep=""))
pheatmap(as.matrix(fpkm_results), legend=TRUE, filename=out_fpkm_heatmap)

################## boxplot of fPKM #############################################
fpkm_data_boxplot <- merge(x=fpkm_data, y=group_table2)
out_fpkm_violin <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_fpkm_violin.pdf", sep=""))
fpkm_violin <- ggplot(data = fpkm_data_boxplot, mapping = aes(x=sample, y=fpkm)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
fpkm_violin_group <- ggplot(data = fpkm_data_boxplot, mapping = aes(x=Group, y=fpkm)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
pdf(out_fpkm_violin)
fpkm_violin
fpkm_violin_group
dev.off()




if (feature_name %in% c("CDS", "gene", "transcript", "exon")){
    
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

    #PCA plot
    out_pca <- file.path(out_dir, paste(strsplit(basename(reads_file), ".tsv")[[1]], "_PCA.pdf", sep=""))
    pdf(out_pca)
    dds_vts <- DESeq2::varianceStabilizingTransformation(dds)
    DESeq2::plotPCA(dds_vts, intgroup = c("Group"))
    dev.off()

    for (n in 1:length(all_pairs) ) {
        # filename strings for each comparisons
        filename <- paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "et.csv", sep = "__")
        filename_sig <- paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "sig.csv", sep = "__")

        # sample matrix
        pair1 <- as.character(all_pairs[[n]][1])
        pair2 <- as.character(all_pairs[[n]][2])
        pairs <- c(pair1, pair2)

        # exact test
        deseq_diff <- DESeq2::results(dds, contrast = c("Group", pair1, pair2))
        deseq_diff <- deseq_diff[order(as.numeric(deseq_diff$pvalue)),]
        deseq_sig <- subset(deseq_diff, pvalue < as.numeric(pcutoff))


        #plot
        out_ma <- file.path(out_dir, paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "MA.pdf", sep = "__"))
        pdf(out_ma)
        plotMA(deseq_diff)
        dev.off()

        # full path to output directory
        out_table <- file.path(out_dir, filename)
        out_table_sig <- file.path(out_dir, filename_sig)

        # write the file
        write.csv(deseq_diff, out_table)
        write.csv(deseq_sig, out_table_sig)

    } 
    } 
