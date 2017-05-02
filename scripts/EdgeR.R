#!/usr/bin/env Rscript
#
#
# author = "Migun Shakya"
# email = "microbeatic@gmail.com"
#
#

# usage: Rscript ../../bin/rpkm_from_reads.R -r sense_reads.csv -g IV_A_PR8_34_mofied_cds.gff -s sum_exp_stats.txt -e ../2016-10-14/ExpDescrNonPolyA_nogz.txt -t region -n segment
#

#
# To calculate the rplm
#


library(optparse)
library(edgeR)
library(rtracklayer)
library(dplyr)


option_list <- list(
  make_option(c("-r", "--reads_table"), action = "store",
              help = "reads table generated from mapped_summary.R"),
  make_option(c("-g", "--gff_file"), action = "store",
              help = "gff file to extract information on length"),
  make_option(c("-s", "--summ_table"), action = "store",
              help = "sum_exp_stats.txt table that contains information on reads mapped"),
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

opt <- parse_args(OptionParser(option_list = option_list))

reads_file <- opt$reads_table
gff <- opt$gff_file
summ <- opt$summ_table
group_file <- opt$exp_desn
feat_name <- opt$name
tot_reads <-opt$tot_reads
feat_type <- opt$type
out_table <- opt$out_table

# read in the reads table
reads_table <- read.delim(reads_file, check.names=FALSE, stringsAsFactors=FALSE, sep=",")
print(reads_table)
print(reads_table[["Sample"]])
# read in the gff table
gff_table <- readGFF(gff)
# print(gff_table)
gff_table <- dplyr::filter(as.data.frame(gff_table), type == feat_type)
gff_table <- dplyr::select(gff_table, get(feat_name), start, end, strand)
print(gff_table)
# gff_table <- subset(gff_table, get(feat_name) %in% reads_table[["Sample"]])
gff_table <- gff_table[gff_table[[feat_name]] %in% reads_table[["Sample"]], ]


# gff_table <- subset(gff_table, protein_id %in% c("NP_040978.1", "NP_040979.2", "NP_040980.1","NP_040981.1", "NP_040982.1", "NP_040983.1"))

# gff_table <- subset(gff_table, protein_id %in% c("NP_040978.1", "NP_040979.2", "NP_040980.1","NP_040981.1", "NP_040982.1", "NP_040983.1"))
#  "NP_040984.1", "NP_040985.1", "NP_040986.1", "NP_040987.1", "YP_006495785.1", "YP_418248.1"))
# gff_table <- dplyr::filter(gff_table, get(feat_name, envir=as.environment(gff_table)) %in% as.character(reads_table$Sample))
# print(gff_table)
gff_table$length <- abs(gff_table$end - gff_table$start) + 1

# TODO: need a way to input variable here (feat_name instead of protein_id)

gff_table <- gff_table %>% 
  group_by(protein_id) %>% 
  summarise(length = sum(length))
print(gff_table)


# remove rows from reads table that are not in gff table
print(feat_type)
print(gff_table[[feat_type]])
reads_table <- dplyr::filter(reads_table,  Sample %in% gff_table[[feat_name]])
print(reads_table)
row.names(reads_table) <- reads_table$Sample
reads_table$Sample <- NULL
print(reads_table)



# read in the summary table with information in reads mapped
summ_table <- read.delim(summ, row.names = NULL)
print(head(summ_table))

# colnames(summ_table) <- c("Sample_ID", "Total_reads", "Total_HighQuality_reads", 
#                          "Total_unmapped_reads", "Total_prokarya_reads",
#                          "unknown1", "Total_prokaryote_reads",
#                          "unknown2","Proper_paired_prokaryote_reads")

summ_table <- select(summ_table, Sample_ID, get(tot_reads))
print(summ_table)

# read in the table with group info
group_table <- read.delim(group_file)
group_table <- select(group_table, ID, group)

print("migun")
DGE_object <- DGEList(counts = reads_table, genes = gff_table, lib.size = summ_table[,2],
                      group = group_table$group)
print("migun")
print(DGE_object)

rpkm_results <- rpkm(DGE_object)
cpm_results <- cpm(DGE_object)


write.csv(rpkm_results, file = paste0(reads_file, "_rpkm.csv"))
write.csv(cpm_results, file = paste0(reads_file, "_cpm.csv"))
