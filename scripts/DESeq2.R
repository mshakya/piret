#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(GenomicRanges)
library(RPiReT)

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


    # Create pairwise comparisons to find DGEs
pair.comb <- function(exp_des){
        # get all pariwise combination from experimental design file
        exp_desn <- read.table(exp_des, sep = "\t", header = TRUE, comment.char = "")
        categories <- unique(exp_desn$Group)
        pairs <- combn(categories, 2, simplify = FALSE)
        return(pairs)
    }

if (feature_name %in% c("CDS", "gene", "transcript", "exon")){
    deseq_ds <- RPiReT::feat2deseq2(feat_count = reads_file, exp_desn = group_file)
    fpkm_table <- RPiReT::DESeq2FPKM(deseq_ds, reads_file, out_dir )
    fpm_table <- RPiReT::DESeq2FPM(deseq_ds, reads_file, out_dir )
    out_fpkm <- file.path(out_dir, paste(strsplit(basename(reads_file),
                                              ".tsv")[[1]],
                                    "_FPKM.csv", sep=""))
    out_fpm <- file.path(out_dir, paste(strsplit(basename(reads_file),
                                              ".tsv")[[1]],
                                    "_FPM.csv", sep=""))
    RPiReT::FPKM_heatmap(out_fpkm, group_table, out_dir)
    RPiReT::FPM_heatmap(out_fpm, group_table, out_dir)
    RPiReT::DESeq2_histogram(out_fpm, group_file, out_dir, "FPM", feature_name)
    RPiReT::DESeq2_histogram(out_fpkm, group_file, out_dir, "FPKM", feature_name)
    RPiReT::DESeq2_violin(out_fpm, group_file, out_dir, "FPM", feature_name)
    RPiReT::DESeq2_violin(out_fpkm, group_file, out_dir, "FPKM", feature_name)
    RPiReT::DESeq2_CAplot(feat_count = reads_file, DESeq2_object = deseq_ds,
                          outdir = out_dir, feature_name = feature_name)
    # calculate size factors
    dds <- DESeq2::DESeq(deseq_ds)
    # get all possible pairs
    all_pairs <- pair.comb(group_file)

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
        summ <- DESeq2::summary.DESeqResults(deseq_diff)
        write.csv(summ, "test_summ.txt")
        deseq_diff <- deseq_diff[order(as.numeric(deseq_diff$pvalue)),]
        deseq_sig <- subset(deseq_diff, padj < as.numeric(pcutoff))

        #plot
        out_ma_pdf <- file.path(out_dir, paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "MA.pdf", sep = "__"))
        pdf(out_ma_pdf)
        plotMA(deseq_diff)
        dev.off()
        out_ma_png <- file.path(out_dir, paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "MA.png", sep = "__"))
        png(out_ma_png)
        plotMA(deseq_diff)
        dev.off()

        RPiReT::DESeq2_summary(object = deseq_diff, alpha = as.numeric(pcutoff),
                    pair1 = pair1, pair2 = pair2, feature_name = feature_name,
                    outdir = out_dir)

        # full path to output directory
        out_table <- file.path(out_dir, filename)
        out_table_sig <- file.path(out_dir, filename_sig)

        # write the file
        write.csv(deseq_diff, out_table)
        write.csv(deseq_sig, out_table_sig)

    } 
}

#TODO: Need to make sure that all gff lines have IDs