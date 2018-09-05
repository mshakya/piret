#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
  make_option(c("-d", "--dir"), action = "store",
              help = "directory containing results of EdgeR or DESeq2"),
  make_option(c("-m", "--method"), action = "store",
              help = "output of which results to conduct gage analysis on"),
  make_option(c("-o", "--out_dir"), action = "store",
              help = "output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

#
dir <- opt$dir
method <- opt$method

#==============================================================================#
# create the output directory
# ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))
#==============================================================================#

# grab fold change files
if ( method == "edgeR"){
    temp = list.files(path=dir, pattern="gene_*_et.csv")

} else if (method == "DESeq2")
{   
    temp = list.files(pattern="gene_*_et.csv")
    print(temp)

}
# read in all the 