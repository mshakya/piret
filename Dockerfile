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
RUN conda install -c bioconda hisat2=2.0.5
RUN conda install -c bioconda faqcs
RUN conda install -c bioconda samtools=1.6
RUN conda install -c bioconda gffread=0.9.12
RUN conda install -c bioconda bamtools=2.4.0
RUN conda install -c bioconda subread=1.6.0
RUN conda install -c bioconda stringtie=1.3.3
