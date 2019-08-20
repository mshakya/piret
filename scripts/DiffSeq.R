#!/usr/bin/env Rscript
#
# usage: Rscript ../../bin/DiffSeq.R -r sense_reads.csv -g IV_A_PR8_34_mofied_cds.gff -s sum_exp_stats.txt -e ../2016-10-14/ExpDescrNonPolyA_nogz.txt -t region -n segment
#


library(optparse)
library(edgeR)
library(pheatmap)
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
group_table <- read.delim(group_file, comment.char = "")
group_table <- select(group_table, X.SampleID, Group)

# max for column position of samples
samp_cols <- as.numeric(nrow(group_table)) + 6

# read count table
reads_table <- read.table(reads_file,
  check.names = FALSE, stringsAsFactors = FALSE, sep = "\t",
  header = TRUE)

# change the column name to just have sample names
col_list <- colnames(reads_table)
colnames(reads_table) <- sapply(col_list, get_samp_name)

# convert DGE object
dge_object <- DGEList(counts = reads_table[, 7:samp_cols],
  genes = reads_table[, 1:6],
  group = group_table$Group)

# remove features that have less than 1 read per million in 3 samples
dge_object <- dge_object[rowSums(1e+06 * dge_object$counts / expandAsMatrix(dge_object$samples$lib.size, dim(dge_object)) > 1) >= 3, ]
# calculate normalization factor
dge_object <- calcNormFactors( dge_object )

# calculate rpkm and cpm
rpkm_results <- rpkm(dge_object)
cpm_results <- cpm(dge_object)


# et <- exactTest(dge_object, pair=levels(group_table$Group))
# topTags(et)
# top <- topTags(et, n=nrow(dge_object$counts))$table
# head(top)




# y <- exactTest(DGE_object)
# print(y)

edge_rpkm_out <- paste0(strsplit(reads_file, ".", fixed=TRUE)[[1]][1], ".edge_rpkm")
write.csv(rpkm_results, file = paste0(strsplit(reads_file, ".", fixed=TRUE)[[1]][1], ".edge_rpkm"))

#RPKM heatmap
edge_rpm_mat <- read.table(edge_rpkm_out, sep = ",", header=TRUE, row.names=1)
print(edge_rpm_mat)


#CPM heatmap
# construct_heatmap()
write.csv(cpm_results, file = paste0(strsplit(reads_file, ".", fixed=TRUE)[[1]][1], ".edge_cpm"))

pdf("~/Desktop/test1.pdf")
plotMDS( dge_object,  main = "MDS Plot for Count Data", labels = colnames( dge_object$counts ) )
plotSmear(dge_object)
edgeR::plotMD.DGEList(dge_object)
dev.off()

rownames(group_table) <- group_table$X.SampleID
group_table$X.SampleID <- NULL
print(group_table)
heat_obj <- pheatmap::pheatmap(edge_rpm_mat, filename = "~/Desktop/test.pdf",
  annotation_col = group_table)