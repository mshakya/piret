# This is the Dockefile to build PiReT (mshakya/PiReT)
# Base Docker Image
FROM continuumio/miniconda3

# Maintainer
MAINTAINER Migun Shakya, migun@lanl.gov

# Update the system
RUN apt-get -y update
RUN apt-get -y install build-essential git-all wget
RUN apt-get clean

# install all piret dependencies
# RUN git clone https://github.com/mshakya/piret.git
# RUN cd piret
RUN conda create -n piret_env python=3.6.6 --yes
RUN conda install -c bioconda faqcs -n piret_env --yes
RUN conda install -c bioconda star hisat2 subread -n piret_env --yes
RUN conda install -c bioconda subread stringtie -n piret_env --yes
RUN conda install -c bioconda samtools bamtools bedtools -n piret_env --yes
RUN conda install -c bioconda diamond=0.9.24 -n piret_env --yes
RUN conda install -c bioconda perl-lwp-protocol-https -n piret_env --yes
RUN conda install -c bioconda perl-json -n piret_env --yes
# RUN source ~/.bashrc
# RUN source activate piret_env
RUN cd thirdparty
RUN rm -rf eggnog-mapper
RUN git clone https://github.com/mshakya/eggnog-mapper.git
# RUN cd eggnog-mapper
# # RUN python download_eggnog_data.py -y
RUN cd ..
# RUN cd ..
RUN Rscript --no-init-file -e "if('BiocManager' %in% rownames(installed.packages()) == FALSE){install.packages('BiocManager',repos='https://cran.r-project.org')}";
RUN Rscript --no-init-file -e "if('optparse' %in% rownames(installed.packages()) == TRUE){packageVersion('optparse');}"  | awk '{print " - found optparse "$2}'
RUN Rscript --no-init-file -e "if('optparse' %in% rownames(installed.packages()) == FALSE){install.packages('optparse',repos='https://cran.r-project.org')}";
RUN Rscript --no-init-file -e "if('tidyverse' %in% rownames(installed.packages()) == TRUE){packageVersion('tidyverse');}"  | awk '{print " - found tidyverse "$2}'
RUN Rscript --no-init-file -e "if('tidyverse' %in% rownames(installed.packages()) == FALSE){install.packages('tidyverse',repos='https://cran.r-project.org')}";
RUN Rscript --no-init-file -e "if('reshape2' %in% rownames(installed.packages()) == TRUE){packageVersion('reshape2');}"  | awk '{print " - found reshape2 "$2}'
RUN Rscript --no-init-file -e "if('reshape2' %in% rownames(installed.packages()) == FALSE){install.packages('reshape2',repos='https://cran.r-project.org')}";
RUN Rscript --no-init-file -e "if('pheatmap' %in% rownames(installed.packages()) == TRUE){packageVersion('pheatmap');}"  | awk '{print " - found pheatmap "$2}'
RUN Rscript --no-init-file -e "if('pheatmap' %in% rownames(installed.packages()) == FALSE){install.packages('pheatmap',repos='https://cran.r-project.org')}";
RUN Rscript --no-init-file -e "if('edgeR' %in% rownames(installed.packages()) == TRUE){packageVersion('edgeR');}"  | awk '{print " - found edgeR "$2}'
RUN Rscript --no-init-file -e "if('edgeR' %in% rownames(installed.packages()) == FALSE){BiocManager::install('edgeR')}";
RUN Rscript --no-init-file -e "if('DESeq2' %in% rownames(installed.packages()) == TRUE){packageVersion('DESeq2');}"  | awk '{print " - found DESeq2 "$2}'
RUN Rscript --no-init-file -e "if('DESeq2' %in% rownames(installed.packages()) == FALSE){BiocManager::install('DESeq2')}";
RUN Rscript --no-init-file -e "if('pathview' %in% rownames(installed.packages()) == TRUE){packageVersion('pathview');}"  | awk '{print " - found pathview "$2}'
RUN Rscript --no-init-file -e "if('pathview' %in% rownames(installed.packages()) == FALSE){BiocManager::install('pathview')}";
RUN Rscript --no-init-file -e "if('gage' %in% rownames(installed.packages()) == TRUE){packageVersion('gage');}"  | awk '{print " - found gage "$2}'
RUN Rscript --no-init-file -e "if('gage' %in% rownames(installed.packages()) == FALSE){BiocManager::install('gage')}";
RUN Rscript --no-init-file -e "if('ballgown' %in% rownames(installed.packages()) == TRUE){packageVersion('ballgown');}"  | awk '{print " - found ballgown "$2}'
RUN Rscript --no-init-file -e "if('ballgown' %in% rownames(installed.packages()) == FALSE){BiocManager::install('ballgown')}";
python setup.py install
