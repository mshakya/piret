library(ballgown)
library(optparse)
#compare the gene expression of control over 

option_list <- list(
  make_option(c("-i", "--bg_folder"), action = "store",
              help = "ballgown folder"),
  make_option(c("-e", "--exp_desn"), action = "store",
              help = "experimental design file"),
    make_option(c("-q", "--q_value"), action = "store",
              help = "q_value to filter differential genes")
)

opt <- parse_args(OptionParser(option_list = option_list))

create.bg <- function(bg_folder, exp_des){
    # returns a bg class object gived bg demanded folder and experimental
    # design file

  dir_names <- list.dirs(path = bg_folder, full.names = TRUE, recursive = FALSE)
    bg <- ballgown::ballgown(samples = dir_names)
    exp_tbl <- read.table(exp_des, sep = "\t", header = TRUE)
    exp_tbl <- subset(exp_tbl, select = c("ID", "Group"))
    colnames(exp_tbl) <- c("id", "group")
    pData(bg) <- exp_tbl
    return(bg)
}

pair.comb <- function(exp_des){
    # get all pariwise combination from experimental design file
    exp_desn <- read.table(exp_des, sep = "\t", header = TRUE)
    categories <- unique(exp_desn$Group)
    pairs <- combn(categories, 2, simplify = FALSE)
    pairs_reg <- as.character(pairs_)
    # pairs_reg <- paste(pairs, sep = "|")
    return(pairs_reg)
}

groups <- pair.comb(opt$exp_desn)
dir_names <- list.dirs(path = opt$bg_folder, 
    full.names = TRUE, recursive = FALSE)

for (grp in groups){
    for (g in as.character(grp)) {
        
        <- grep(g, dir_names, value = TRUE)

    exp_desn <- read.table(exp_des, sep="\t", header=TRUE)
    subset()
    dir_names <- list.dirs(path = bg_folder, full.names = TRUE, recursive = FALSE)
    
    bg <- ballgown::ballgown(samples = dir_names)

    }

    dir_names <- list.dirs(path = opt$bg_folder,
        full.names = TRUE, recursive = FALSE)
    sampes = 

}

comp_pair <- function(bg_){
    # conduct all possible stattest pairwise comparison 

    dir_names <- list.dirs(path = bg_folder, full.names = TRUE, recursive = FALSE)
    bg <- ballgown::ballgown(samples = )

}


bg <- create.bg(bg_folder = opt$bg_folder, exp_des = opt$exp_desn)

diff_genes <- stattest(bg, feature = "gene",
                       meas = "FPKM", covariate = "group")

sig_transcripts <- subset(diff_genes, qval < opt$q_value)
sig_transcripts


# bg <- ballgown(samples = c('~/Desktop/RNASeq_Yersinia/2017-04-26/yersinia/ballgown/euk/CTHP1/', '~/Desktop/RNASeq_Yersinia/2017-04-26/yersinia/ballgown/euk/CTHP2', 
#                            '~/Desktop/RNASeq_Yersinia/2017-04-26/yersinia/ballgown/euk/ITHP1', '~/Desktop/RNASeq_Yersinia/2017-04-26/yersinia/ballgown/euk/ITHP2'))
# pData(bg) <- data.frame(id=sampleNames(bg), group=c(0,0,1,1))
# diff_genes <- stattest(bg, feature = 'transcript', meas='FPKM', covariate = 'group')
# subset(diff_genes, qval < 0.001)ok
# #MBNL1
# plotTranscripts('gene10508', bg, 
#                 samples=c('CTHP1', 'CTHP2', 'ITHP1', 'ITHP2'), 
#                 meas='FPKM', colorby='transcript')

# plotMeans('gene10508', bg, groupvar = 'group',
#                 meas='FPKM', colorby='transcript')

# #KIAA1191
# plotTranscripts('gene15641', bg, 
#                 samples=c('CTHP1', 'CTHP2', 'ITHP1', 'ITHP2'), 
#                 meas='FPKM', colorby='transcript')
# #DHX16
# plotTranscripts('gene16680', bg, 
#                 samples=c('CTHP1', 'CTHP2', 'ITHP1', 'ITHP2'), 
#                 meas='FPKM', colorby='transcript')

# #KIF24
# plotTranscripts('gene23667', bg, 
#                 samples=c('CTHP1', 'CTHP2', 'ITHP1', 'ITHP2'), 
#                 meas='FPKM', colorby='transcript')

# #MIR3198-2
# plotTranscripts('gene31307', bg, 
#                 samples=c('CTHP1', 'CTHP2', 'ITHP1', 'ITHP2'), 
#                 meas='FPKM', colorby='transcript')

# ?plotTranscripts
