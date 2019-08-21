#!/usr/bin/env Rscript
library(ballgown)
library(optparse)
#compare the gene expression of control over 

option_list <- list(
    make_option(c("-i", "--bg_folder"), action = "store",
              help = "ballgown folder"),
    make_option(c("-e", "--exp_desn"), action = "store",
              help = "experimental design file"),
    make_option(c("-p", "--p_value"), action = "store",
              help = "p_value to filter differential genes"),
    make_option(c("-n", "--name"), action = "store",
              help = "name of feature from gff file that was chosen to represent each feature"),
    make_option(c("-o", "--outfolder"), action = "store",
              help = "output folder")

)

opt <- parse_args(OptionParser(option_list = option_list))

bg_folder <- opt$bg_folder
pcutoff <- opt$p_value
feature_name <- opt$name
group_file <- opt$exp_desn
out_dir <- opt$outfolder


exp_df <- read.table(group_file, sep = "\t", header = TRUE,
                     comment.char="")

# Functions
create.bg <- function(bg_folder, exp_df){
    # returns a bg class object gived bg demanded folder and experimental
    # design file

    dir_names <- list.dirs(path = bg_folder, full.names = TRUE, recursive = FALSE)
    bg <- ballgown::ballgown(samples = dir_names)
    # exp_df <- subset(exp_df, select = c("ID", "Group"))
    # colnames(exp_df) <- c("id", "group")
    pData(bg) <- exp_df
    return(bg)
}

pair.comb <- function(exp_des){
    # get all pariwise combination from experimental design file
    exp_desn <- read.table(exp_des, sep = "\t", header = TRUE,
                           comment.char="")
    categories <- unique(exp_desn$Group)
    pairs <- combn(categories, 2, simplify = FALSE)
    return(pairs)
}


subset.exp_des <- function(exp_df, pair_combs){
    # split the experimental design file to list of dataframe
    # with all possible combination
    exp_tbl <- subset(exp_df, select = c("X.SampleID", "Group"))
    colnames(exp_tbl) <- c("id", "group")
    i <- 1
    l <- list()
    for (pair in pair_combs){
        pair_exp_df <- subset(exp_tbl, group %in% pair)
        l[[i]] <- pair_exp_df
        i <- i + 1
    }
    return(l)
}



# get groups
groups <- pair.comb(group_file)
list_of_eds <- subset.exp_des(exp_df, groups)
print(list_of_eds)

# list all the folders in ballgown folder
dir_names <- list.dirs(path = opt$bg_folder, 
    full.names = TRUE, recursive = FALSE)

# create the output directory
ifelse(!dir.exists(out_dir), dir.create(out_dir), print("already exist"))

if (feature_name %in% c("gene", "CDS", "transcript")){
    for (df in list_of_eds){
        filename <- paste(unique(df$group)[1], unique(df$group)[2], feature_name, "et.csv", sep = "__")
        filename_sig <- paste(unique(df$group)[1], unique(df$group)[2], feature_name, "sig.csv", sep = "__")
        df$group <- as.factor(df$group)
        bg <- create.bg(bg_folder = bg_folder, exp_df = df)
        print(bg)
        diff_genes <- stattest(bg, feature = feature_name,
                                meas = "FPKM", covariate = "group", getFC=TRUE)
        sig_genes <- subset(diff_genes, pval < pcutoff)
        # full path to output directory
        out_table <- base::file.path(out_dir, filename)
        out_table_sig <- base::file.path(out_dir, filename_sig)
        # write the file
        print(out_table)
        write.csv(diff_genes, out_table)
        write.csv(sig_genes, out_table_sig)
    }
    
}
