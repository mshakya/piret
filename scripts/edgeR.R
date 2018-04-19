#!/usr/bin/env Rscript

library(optparse)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(gridExtra)


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

#
reads_file <- opt$reads_table
pcutoff <- opt$p_cutoff
feature_name <- opt$name
group_file <- opt$exp_desn
out_dir <- opt$out_dir

# read the output of featureCounts
read.counts <- read.table(reads_file, sep = "\t", header = TRUE, comment.char = "#")

# rename the column headers
names(read.counts) <- gsub(".*mapping_results.", "", names(read.counts),
                           perl = TRUE)
names(read.counts) <- gsub("_srt.bam", "", names(read.counts), perl = TRUE)


# convert first column to row names
row.names(read.counts) <- read.counts[, 1]

# get gene information
gene.info <- read.counts[, c(1:6)]

# read the experimental design file
group_table <- read.delim(group_file, row.names = 1)
group_table <- select(group_table, Group)


#TODO: move this piece of code to separate R code so that its done everytime
#rearrange/clean the table and rewrite the column to display in EDGE
outfile <- paste(strsplit(reads_file ,".tsv"), "_sorted.csv", sep="")
read.counts.sort <- read.counts
col_len <- ncol(read.counts.sort)
read.counts.sort$sum <- rowSums(read.counts.sort[7:col_len])
read.counts.sort <- dplyr::arrange(read.counts.sort, desc(sum))
read.counts.sort$sum <- NULL
read.counts.sort <- read.counts.sort[c(base::colnames(gene.info), base::row.names(group_table))]
write.csv(read.counts.sort, file=outfile )


# Reformat summary files for EDGE display
summary_file <- paste(reads_file, ".summary", sep="")
summary_outfile <- paste(strsplit(reads_file ,".summary"), "_summary.csv", sep="")
summary.counts <- read.table(summary_file, sep = "\t", header = TRUE, comment.char = "#")
names(summary.counts) <- gsub(".*mapping_results.", "", names(summary.counts),
                           perl = TRUE)
names(summary.counts) <- gsub("_srt.bam", "", names(summary.counts), perl = TRUE)
summary.counts <- summary.counts[c("Status", base::row.names(group_table))]
write.csv(summary.counts, file=summary_outfile )


# remove first six columns that have metadata
read.counts <- read.counts[, -c(1:6)]
read.counts.non0<- dplyr::filter_all(read.counts, any_vars(. != 0))

# create the output directory
ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))

out_heatmap <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_heatmap.pdf", sep=""))
pdf(out_heatmap)
pheatmap(as.matrix(read.counts.non0), legend=TRUE)
dev.off()

edger_dge <- edgeR::DGEList(counts = read.counts, group = group_table$Group,
                            remove.zeros = FALSE, genes = gene.info)

############### calculate RPKM and CPM #########################################
rpkm_results <- edgeR::rpkm(edger_dge)
cpm_results <- edgeR::cpm(edger_dge)
out_rpkm <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_RPKM.csv", sep=""))
out_cpm <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_CPM.csv", sep=""))
write.csv(rpkm_results, file = out_rpkm)
write.csv(cpm_results, file = out_cpm)
################################################################################

################ histogram of count per million ################################
out_cpm_hist <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_cpm_histograpm.pdf", sep=""))
cpm_results <- dplyr::filter_all(as.data.frame(cpm_results), any_vars(. != 0))
cpm_data <- reshape2::melt(as.data.frame(cpm_results), variable.name="sample", value.name="CPM")
cpm_hist <- ggplot(data = cpm_data, mapping = aes(x = CPM)) +  theme_bw() +
        geom_histogram(bins=100) + xlab("CPM") + ylab(feature_name) + facet_wrap(~sample)
ggsave(out_cpm_hist, cpm_hist, device = "pdf")

################## boxplot of count per million ################################
group_table2 <- add_rownames(group_table, "sample")
cpm_data_boxplot <- merge(x=cpm_data, y=group_table2)
out_cpm_violin <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_cpm_violin.pdf", sep=""))
cpm_violin <- ggplot(data = cpm_data_boxplot, mapping = aes(x=sample, y=CPM)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
cpm_violin_group <- ggplot(data = cpm_data_boxplot, mapping = aes(x=Group, y=CPM)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
pdf(out_cpm_violin)
cpm_violin
cpm_violin_group
dev.off()

##############################histogram of rpkm ################################
out_rpkm_hist <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_rpkm_histograpm.pdf", sep=""))
rpkm_results <- dplyr::filter_all(as.data.frame(rpkm_results), any_vars(. != 0))
rpkm_data <- reshape2::melt(as.data.frame(rpkm_results), variable.name="sample", value.name="rpkm")
rpkm_hist <- ggplot(data=rpkm_data, mapping=aes(x=rpkm)) +  theme_bw() +
        geom_histogram(bins=100) + xlab("rpkm") + ylab(feature_name) + facet_wrap(~sample)
ggsave(out_rpkm_hist, rpkm_hist, device = "pdf")

##############################heatmap of rpkm ################################
out_rpkm_heatmap <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_rpkm_heatmap.pdf", sep=""))
pheatmap(as.matrix(rpkm_results), legend=TRUE, filename=out_rpkm_heatmap)

################## boxplot of RPKM #############################################
rpkm_data_boxplot <- merge(x=rpkm_data, y=group_table2)
out_rpkm_violin <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_rpkm_violin.pdf", sep=""))
rpkm_violin <- ggplot(data = rpkm_data_boxplot, mapping = aes(x=sample, y=rpkm)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
rpkm_violin_group <- ggplot(data = rpkm_data_boxplot, mapping = aes(x=Group, y=rpkm)) +  theme_bw() +
        geom_violin(aes(fill = factor(Group)))
pdf(out_rpkm_violin)
rpkm_violin
rpkm_violin_group
dev.off()

if (feature_name %in% c("CDS", "gene", "transcript", "exon")){

    if (0 %in% colSums(edger_dge$counts)) {
        print("One of the sample does not have any reads mapped to it, so no further analysis will be done!")
    } else {
        # only keep genes/features that have an average of 1 CPM
        keep <- rowSums(edgeR::cpm(edger_dge)) >= nrow(group_table)

        edger_dge <- edger_dge[keep, ]
        if (nrow(edger_dge$counts) < 10 ) {
            print("Not enough data to calculate dispersion factors!")
        } else {
        # recompute library size after filtering
        edger_dge$samples$lib.size <- colSums(edger_dge$counts)
        # calculate the normalization factors for the library size
        edger_dge <- edgeR::calcNormFactors(edger_dge, method = "TMM")
        # estimate the dispersion for all read counts across all samples
        edger_dge <- edgeR::estimateDisp(edger_dge, design.robust = TRUE)

        # plot MDS
        out_mds_hist <- file.path(out_dir, paste(strsplit(basename(reads_file), ".csv")[[1]], "_MDS.pdf", sep=""))
        pdf(out_mds_hist)
        limma::plotMDS(edger_dge, top = 1000, labels = edger_dge$sample$ID,
               main = "edgeR MDS Plot")
        dev.off()

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
            filename <- paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "et.csv", sep = "__")
            filename_sig <- paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name,"sig.csv", sep = "__")

            # sample matrix
            pair1 <- as.character(all_pairs[[n]][1])
            pair2 <- as.character(all_pairs[[n]][2])
            pairs <- c(pair1, pair2)

          # exact test
            edger_et <- edgeR::exactTest(edger_dge, pair = pairs)
            edger_et_table <- edgeR::topTags(edger_et, n = nrow( edger_et$table ), sort.by = "PValue" )$table
            edger_et_sigtable <-  subset(edger_et_table, PValue < pcutoff)
           
           # write summary table
            out_file_summ <- file.path(out_dir, paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "summary.csv", sep = "__"))
            summ_table <- summary(decideTests(edger_et))
            colnames(summ_table) <- paste(colnames(summ_table), feature_name, sep="_")
            write.csv(summ_table, out_file_summ)
            
            #plot
            out_md <- file.path(out_dir, paste(all_pairs[[n]][1], all_pairs[[n]][2], feature_name, "MD.pdf", sep = "__"))
            pdf(out_md)
            plotMD(edger_et)
            dev.off()
          
          # full path to output directory
            out_table <- file.path(out_dir, filename)
            out_table_sig <- file.path(out_dir, filename_sig)

          # write the file
            write.csv(edger_et_table, out_table)
            write.csv(edger_et_sigtable, out_table_sig)
            }

    }
}
}
