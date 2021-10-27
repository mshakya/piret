#!/usr/bin/env bash

if [ "$1" != "" ]; then
    echo "Environment name will be $1"
    env=$1
else
    echo "Positional parameter 1 is empty"
    env="piret_env"
fi


set -e # Exit as soon as any line in the bash script fails

ROOTDIR=$( cd $(dirname $0) ; pwd -P ) # path to main PiReT directory


echo
exec &> >(tee -a  install.log)
exec 2>&1 # copies stderr onto stdout

conda create -n $env python=3.6.6 --yes
conda install -c bioconda faqcs -n $env --yes
conda install -c bioconda star hisat2 subread -n $env --yes
conda install -c bioconda subread stringtie -n $env --yes
conda install -c bioconda samtools bamtools bedtools -n $env --yes
conda install -c bioconda diamond=0.9.24 -n $env --yes
conda install -c bioconda perl-lwp-protocol-https -n $env --yes
conda install -c bioconda perl-json -n $env --yes
conda install -c conda-forge time -n $env --yes
conda install -c r r-tidyverse r-reshape2 r-optparse r-pheatmap
conda install -c bioconda bioconductor-deseq2 bioconductor-edger bioconductor-gage bioconductor-ballgown
export PATH=/opt/conda/envs/piret_env/bin:$PATH
source activate $env
rm -rf thirdparty/eggnog-mapper
git clone https://github.com/mshakya/eggnog-mapper.git thirdparty/eggnog-mapper
# R packages
echo "if('BiocManager' %in% rownames(installed.packages()) == FALSE) {install.packages('BiocManager', repos='https://cran.r-project.org')}" | Rscript -
#echo "if('optparse' %in% rownames(installed.packages()) == TRUE) {packageVersion('optparse')}"  | Rscript -
#echo "if('optparse' %in% rownames(installed.packages()) == FALSE) {install.packages('optparse',repos='https://cran.r-project.org')}"  | Rscript -
# install tidyverse
#echo "if('tidyverse' %in% rownames(installed.packages()) == TRUE){packageVersion('tidyverse')}"  | Rscript -
#echo "if('tidyverse' %in% rownames(installed.packages()) == FALSE){install.packages('tidyverse',repos='https://cran.r-project.org')}"  | Rscript -
# install R reshape2 packages
# echo "if('reshape2' %in% rownames(installed.packages()) == TRUE){packageVersion('reshape2')}"  | Rscript -
# echo "if('reshape2' %in% rownames(installed.packages()) == FALSE){install.packages('reshape2',repos='https://cran.r-project.org')}"  | Rscript -
# install R pheatmap packages
#echo "if('pheatmap' %in% rownames(installed.packages()) == TRUE){packageVersion('pheatmap')}"  | Rscript -
#echo "if('pheatmap' %in% rownames(installed.packages()) == FALSE){install.packages('pheatmap',repos='https://cran.r-project.org')}"  | Rscript -
# install R edgeR packages
#echo "if('edgeR' %in% rownames(installed.packages()) == TRUE){packageVersion('edgeR')}"  | Rscript -
#echo "if('edgeR' %in% rownames(installed.packages()) == FALSE){BiocManager::install('edgeR')}"  | Rscript -
# install R deseq2 packages
#echo "if('DESeq2' %in% rownames(installed.packages()) == TRUE){packageVersion('DESeq2')}"  | Rscript -
#echo "if('DESeq2' %in% rownames(installed.packages()) == FALSE){BiocManager::install('DESeq2')}"  | Rscript -
# install R gage package
#echo "if('gage' %in% rownames(installed.packages()) == TRUE){packageVersion('gage')}"  | Rscript -
#echo "if('gage' %in% rownames(installed.packages()) == FALSE){BiocManager::install('gage')}"  | Rscript -
# install R ballgown package
#echo "if('ballgown' %in% rownames(installed.packages()) == TRUE){packageVersion('ballgown')}"  | Rscript -
#echo "if('ballgown' %in% rownames(installed.packages()) == FALSE){BiocManager::install('ballgown')}"  | Rscript -
python setup.py install
pytest --cov=piret tests/



echo "
All done!
Run
source activate $env
piret -h
for usage.
Read the README for more information!
	"
