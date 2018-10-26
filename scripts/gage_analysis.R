#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
  make_option(c("-d", "--dir"), action = "store",
              help = "directory containing results of EdgeR or DESeq2"),
  make_option(c("-m", "--method"), action = "store",
              help = "output of which results to conduct gage analysis on"),
  make_option(c("-c", "--org_code"), action = "store",
              help = "kegg orgainsm code to download pathway information"),
  make_option(c("-o", "--out_dir"), action = "store",
              help = "output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

#
dir <- opt$dir
method <- opt$method
code <- opt$org_code

#==============================================================================#
# create the output directory
# ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))
#==============================================================================#

# grab fold change files
if ( method == "edgeR"){
    fcs = list.files(path = dir, pattern = "gene_*_et.csv", full.names = TRUE)

} else if (method == "DESeq2")
{   
    fcs = list.files(path = dir, pattern = "gene_*_et.csv", full.names = TRUE)
    
}


# perform gage analysis, one fc at a time

for (fc in fcs){
    fc_df <- read.csv(fc, row.names = 1)
    print(head(fc_df))

}


gene_sets <- kegg.gsets(species = code)