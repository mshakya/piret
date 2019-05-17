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
RUN conda install -c bioconda piret

# Get latest piret
RUN git clone https://github.com/mshakya/piret.git /root/piret
RUN cd /root/piret
RUN export PATH="/root/piret/bin:$PATH"