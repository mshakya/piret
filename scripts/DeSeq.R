library(Deseq2)

deseq_matrix <- DESeqDataSetFromHTSeqCount(countData = countdata,
    colData = coldata, 

# cmd_args = commandArgs(TRUE)
# pcutoff <- as.numeric(cmd_args[1])
# #mydir <- cmd_args[1]

# targets.ge48h <- read.table(file = "./readcounts.expriment.txt",sep = "\t",  stringsAsFactors = FALSE,  header=TRUE)

# #deseq2
# library(DESeq2)
# sampleTable <- data.frame(sampleName = targets.ge48h[,1], fileName = targets.ge48h$files, condition = targets.ge48h$group)
# dds=DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/", design= ~ condition)
# dds <- DESeq(dds)

# fileout_pca=paste('./figures/Deseq2_PCA', ".pdf", sep="")
# pdf(fileout_pca)
# print(plotPCA(varianceStabilizingTransformation(dds),ntop=1000, intgroup=c("condition")))
# dev.off()

# #EdgeR
# #library(edgeR)

# #eds <- readDGE(targets.ge48h,  comment.char = "!")
# #eds <- estimateCommonDisp(eds, verbose=TRUE)
# #eds <- estimateTagwiseDisp(eds)

# #fileout_bcv=paste('./figures/EdgeR_BCV', ".pdf", sep="")
# #pdf(fileout_bcv)
# #plotBCV(eds, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance")
# #dev.off()


# #fileout_mds=paste('./figures/EdgeR_MDS', ".pdf", sep="")
# #pdf(fileout_mds)
# #plotMDS(eds, top=1000, labels=eds$sample$ID, main="edgeR MDS Plot")
# #dev.off()

# #loop over pairs

# pairs <- read.table(file = "./Deseq_EdgeRpairs.txt",sep = "\t", stringsAsFactors = FALSE,  header=TRUE)
# pair1 <- pairs$group1
# pair2 <-pairs$group2

# for (n in 1:length(pair1) ) {
# if (length(pair1) ==0) next

# filename=paste(pair1[n],pair2[n], sep="__")
# filename_sig=paste(pair1[n],pair2[n],"sig", sep="__")

# #deseq
# fileout_Deseq=paste('./Deseq/Deseq.',filename, ".txt", sep="")
# fileout_sig_Deseq=paste('./Deseq/Deseq.',filename_sig, ".txt", sep="")

# filename_deseq <- results(dds, contrast=c("condition",pair1[n],pair2[n]))
# if (nrow(filename_deseq)==0) next
# filename_deseq <-filename_deseq[order(as.numeric(filename_deseq$pvalue)),]
# write.table(filename_deseq,sep = "\t",  file=fileout_Deseq)

# #significant
# if (nrow(filename_deseq)==0) next
# filename_sig <-  subset(filename_deseq,filename_deseq$pvalue<pcutoff)
# write.table(filename_sig,sep = "\t",  file=fileout_sig_Deseq)

# #edgeR

# #fileout_EdgeR=paste('./EdgeR/EdgeR.',filename, ".txt", sep="")
# #fileout_sig_EdgeR=paste('./EdgeR/EdgeR.',filename_sig, ".txt", sep="")

# #filename_EdgeR <- exactTest(eds, pair=c(pair1[n],pair2[n]))
# #filename_EdgeR<-topTags( filename_EdgeR, n = nrow( filename_EdgeR$table ) , sort.by = "PValue" )$table
# #if (nrow(filename_EdgeR)==0) next
# #filename_EdgeR <-filename_EdgeR[order(as.numeric(filename_EdgeR$PValue)),]
# #write.table(filename_EdgeR,sep = "\t",  file=fileout_EdgeR)
# #filename_sig_EdgeR <-  subset(filename_EdgeR,filename_EdgeR$PValue<pcutoff)




# #significant
# #if (nrow(filename_sig_EdgeR)==0) next
# #filename_sig_EdgeR <-  subset(filename_EdgeR,filename_EdgeR$PValue<pcutoff)
# #write.table(filename_sig_EdgeR,sep = "\t",  file=fileout_sig_EdgeR)
# #dev.off()
# }
