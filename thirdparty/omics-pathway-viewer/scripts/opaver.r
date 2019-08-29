#!/usr/bin/env Rscript

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
opaver.r --wd <DIR> [OPTIONS]
      --wd=<DIR>     working directory (default: .)

[OPTIONS]
      --pwy=<PATH>   path of exp_pathway.txt
      --exp=<PATH>   path of exp_gene.txt
      --cpd=<PATH>   path of exp_cpd.txt
      --help         print this help
 
Example:
	  ./opaver.r --wd=output\n\n")

	q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

pwyFile <- argsL$pwy
expFile <- argsL$exp
cpdFile <- argsL$cpd
dir     <- argsL$wd

## set outdir as working dir by default
if(is.null(dir)) {
	dir <- getwd()
}
setwd(dir)

## load pathway, gene expression data and metabolic data 
if(is.null(argsL$pwy)) {
	pwyFile <- "exp_pathway.txt"
}
if(is.null(argsL$exp)) {
	expFile <- "exp_gene.txt"
}
if(is.null(argsL$exp)) {
	cpdFile <- "exp_cpd.txt"
}

library(pathview) 
options(bitmapType="cairo")

pathway<-read.table(file=pwyFile, sep="\t", colClasses="character", quote="") 
gene_matrix_table<-read.table(file=expFile, sep="\t", header=TRUE)
cpd_matrix_table<-read.table(file=cpdFile, sep="\t", header=TRUE, quote="")

#only keep sample logFC columns
gene_matrix <- as.matrix(gene_matrix_table[,4:length(gene_matrix_table[1,])])
cpd_matrix  <- as.matrix(cpd_matrix_table[,3:length(cpd_matrix_table[1,])])
n <- length(gene_matrix[1,])-1
cpd_matrix  <- cbind(cpd_matrix, replicate(n, cpd_matrix[,1])) #duplicate cloumn to match the sample size
#gene_matrix[,1] <- log2(gene_matrix[,1])

#proxy settings
#proxy=Sys.getenv("http_proxy")
#Sys.setenv(http_proxy=proxy)

#normalize by 3rd Quantile.
qpp <- quantile(gene_matrix[,1][gene_matrix[,1]>0], na.rm=T, c(0.75))
qpn <- quantile(gene_matrix[,1][gene_matrix[,1]<0], na.rm=T, c(0.25))
qtp <- quantile(gene_matrix[,2][gene_matrix[,2]>0], na.rm=T, c(0.75))
qtn <- quantile(gene_matrix[,2][gene_matrix[,2]<0], na.rm=T, c(0.25))
qcp <- quantile(cpd_matrix[,1][cpd_matrix[,1]>0],   na.rm=T, c(0.75))
qcn <- quantile(cpd_matrix[,1][cpd_matrix[,1]<0],   na.rm=T, c(0.25))


qval <- data.frame(qpp,qpn,qtp,qtn,qcp,qcn)
qval <- abs(qval)
qval[is.na(qval)] <- 1
gene_matrix[,1][which(gene_matrix[,1]>0)] <- gene_matrix[,1][which(gene_matrix[,1]>0)]/max(qval$qpp,qval$qpn)
gene_matrix[,1][which(gene_matrix[,1]<0)] <- gene_matrix[,1][which(gene_matrix[,1]<0)]/max(qval$qpp,qval$qpn)
gene_matrix[,2][which(gene_matrix[,2]>0)] <- gene_matrix[,2][which(gene_matrix[,2]>0)]/max(qval$qtp,qval$qtn)
gene_matrix[,2][which(gene_matrix[,2]<0)] <- gene_matrix[,2][which(gene_matrix[,2]<0)]/max(qval$qtp,qval$qtn)

if( length(cpd_matrix) ){
	cpd_matrix[,1][which(cpd_matrix[,1]>0)] <- cpd_matrix[,1][which(cpd_matrix[,1]>0)]/max(qval$qcp,qval$qcn)
	cpd_matrix[,1][which(cpd_matrix[,1]<0)] <- cpd_matrix[,1][which(cpd_matrix[,1]<0)]/max(qval$qcp,qval$qcn)
	cpd_matrix[,2][which(cpd_matrix[,2]>0)] <- cpd_matrix[,2][which(cpd_matrix[,2]>0)]/max(qval$qcp,qval$qcn)
	cpd_matrix[,2][which(cpd_matrix[,2]<0)] <- cpd_matrix[,2][which(cpd_matrix[,2]<0)]/max(qval$qcp,qval$qcn)
}

#KEDGG ID as rownames
rownames(gene_matrix)=gene_matrix_table[,1]
rownames(cpd_matrix)=cpd_matrix_table[,1]

sample_size<-length(gene_matrix[1,])
for (i in 1:length(pathway$V1)){
	tryCatch(
		pathview(
			pathway.id = pathway$V1[i],
			gene.data = gene_matrix[,1:sample_size],
			cpd.data = cpd_matrix[,1:sample_size],
			gene.idtype = "kegg",
			cpd.idtype = "kegg",
			kegg.native = T,
			species="ko",
			same.layer=T
		),
		error = function(e){
			print( paste("[ERROR] Can't generate map", pathway$V1[i], ":", e) )
			return(NULL)
		}
	)
}
