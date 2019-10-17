<!-- [![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://bioconda.github.io/recipes/viral-ngs/README.html) -->

[![Build Status](https://travis-ci.org/mshakya/piret.svg?branch=master)](https://travis-ci.org/mshakya/piret)
[![codecov](https://codecov.io/gh/mshakya/piret/branch/master/graph/badge.svg)](https://codecov.io/gh/mshakya/piret)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/piret/README.html)
<!-- [![Coverage Status](https://coveralls.io/repos/github/mshakya/piret/badge.svg)](https://coveralls.io/github/mshakya/piret) -->

<!-- [ ![Codeship Status for mshakya/PyPiReT](https://app.codeship.com/projects/56b51ac0-1899-0135-279d-22ebbf49dae4/status?branch=master)](https://app.codeship.com/projects/219002) -->
<!-- [![Code Health](https://landscape.io/github/mshakya/PiReT/master/landscape.svg?style=flat&badge_auth_token=4d9548f2284e4b29b79e1bb813573a87)](https://landscape.io/github/mshakya/PiReT/master) -->


# PiReT

Pipeline for Reference based Transcriptomics.

## 0.0 Installing PiReT

### 0.0.1 Installing directly from conda


<!-- PiReT can be installed using conda:
```
conda install piret
```

Note that this doesnt -->

Coming soon!

### 0.0.2

Installing separatelyt

For getting KO ids for genes, PiReT uses [emapper](https://github.com/eggnogdb/eggnog-mapper). The conda install of PiReT also includes emapper. However, its database need to be downloaed following instruction [here](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2). Briefly,

<!-- PiReT uses bioinformatic tools, many of which are available in [bioconda](https://bioconda.github.io). For installing `PiReT` we have provided a script `INSTALL.sh` that checks for required dependencies (including their versions) are installed and in your path, and installs it in directories within `PiReT` if not found. Additionally, `sudo` privileges are not needed for installation. A log of all installation can be found in `install.log` -->

### 0.1 Dependencies
PiReT requires following dependencies, all of which should be installed and in the PATH.

#### 0.1.0 Programming/Scripting languages
- [Python >=v3.6.3](https://www.python.org/downloads/release/python-2712/)
    - The pipeline is not compatible with Python v3.0 or higher.
- [R >=v3.3.1](https://www.r-project.org)
- [Perl >=v5.26.2](https://www.perl.org/)

#### 0.1.1 Installing dependencies
- [conda v4.2.13](http://conda.pydata.org/docs/index.html)
    If conda is not installed, `INSTALL.sh` will download and install [miniconda](http://conda.pydata.org/miniconda.html), a "mini" version of `conda` that only installs handful of packages compared to [anaconda](https://docs.continuum.io/anaconda/pkg-docs)


#### 0.1.2 Third party softwares/packages
- [samtools (>=v1.6)](http://www.htslib.org)
- [HiSat2 (>=v2.1.0)](https://ccb.jhu.edu/software/hisat/index.shtml)
- [featurecount (>=v1.6.3)](https://academic.oup.com/bioinformatics/article/30/7/923/232889/featureCounts-an-efficient-general-purpose-program)
- [stringTie (>=v1.3.4d)](https://ccb.jhu.edu/software/stringtie/)

#### 0.1.3 R packages
- [edgeR (>=v3.14.0)](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [DEseq2 (>=v1.12.4)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [ballgown (>=v2.8.0)](https://bioconductor.org/packages/release/bioc/html/ballgown.html)

#### 0.1.4 Python packages
- [luigi (>=v2.6.1)](https://github.com/spotify/luigi)
- [pandas (>=v0.19.2)](http://pandas.pydata.org/)
- [plumbum (>=v1.6.3)](https://plumbum.readthedocs.io/en/latest/)
- [Biopython (>=v1.68)](https://github.com/biopython/biopython.github.io/)
- [gffread (>=v0.8.4rc1)](https://pythonhosted.org/gffutils/)


## Installtion

```
  conda create -n piret_env python=3.6.6
  conda install -c bioconda faqcs star samtools bamtools hisat2 subread -n piret_env
  conda install -c r r -n piret_env

```
## 1.0 Test
We have provided test data set to check if the installation was successful or not. `fastq` files can be found in `tests/fastqs` and corresponding reference fasta files are found in `tests/data`. To run the test, from within `piret` directory:


For running tests on eukaryote datasets:

```
$LUIGI_CONFIG_PATH="full_path_to/piret/tests/test_euk.cfg" bin/piret -c tests/test_euk.cfg -d tests/test_euk -e tests/test_euk.txt
```

For running tests on prokarya datasets:

```
$LUIGI_CONFIG_PATH="full_path_to/piret/tests/test_prok.cfg" bin/piret -c tests/test_prok.cfg -d tests/test_prok -e tests/test_prok.txt
```

For running tests using `both` prokarya and eukarya datasets:
```
$LUIGI_CONFIG_PATH="full_path_to/piret/tests/test_both.cfg" bin/piret -c tests/test_prok.cfg -d tests/test_prok -e tests/test_both.txt
```


## 2.0 Running PiReT
```
usage: piret [-h] -d WORKDIR -e EXPDSN -c CONFIG [-v]

piret

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

required arguments:
  -d WORKDIR            working directory where all output files will be
                        processed and written (default: None)
  -e EXPDSN             tab delimited experimental design file
  -c CONFIG, --config CONFIG
                        luigi config file for setting parameters that control
                        each step, see github repo for an example (default:
                        None)

Example runs:

        piret -d <workdir> -e <design file>  -c <config file>
```

### 2.1 Experimental design file
An experimental design file consist of sample name (SampleID), full path to fastq files (Files), and different groups of your samples (Group). We recommend that you use a text editor like BBedit or TextWrangler to generate the tab delimited experimental design file. Exporting a tab delimited file directly from Excel tend to cause formatting problem. If possible, please avoid any special characters in sample names and group names.
  
  For example:
  ```
  samp1, samp_1 : good name
  samp 1, samp.1: not a good name and will likely cause errors.

  ```
  A sample of experimental design file can be found [here](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/experimental_design.txt). 

### 2.2 Config file
All options are set in the config file.


## 3.0 OUTPUT

All the outputs will be within the `working directory`. The main output file is a concatenated JSON file called `out.json`.



- `samp2`: The name of this directory corresponds to sample name. Within this folder there are two sub-folders:

  - `mapping_results`
    This folder contains reads mapped using *hisat2* in following formats. If `splice_sites_gff.txt` is present, **hisat2** aligns based on known splice sites.
    - `*.sam`: outputs of *hisat2*
    - `*.bam`: generated from `.sam`
    - [mapping.log](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/mapping.log): Alignment summary file from `hisat2`.
    - [`*sTie.tab`](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/samp3_prok_sTie.tab): Tab delimited file with Coverage, FPKM, TPM, for all the genes and novel transcripts. Generated using string tie.
    - [`*sTie.gtf`](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/samp3_prok_sTie.gtf): Primay GTF formatted output of stringtie.
  - `trimming_results`
  This folder contains results of quality trimming and filtering using FaQC.
    - [`*_qc_report.pdf`](https://github.com/mshakya/PyPiReT/blob/master/examples/samp3_qc_report.pdf): A QC report file with figures.
    - [fastqCount.txt](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/fastqCount.txt): A text file with summary of read counts.
    - *trimmed.fastq: Pair of trimmed fastq files.
    - *unpaired.trimmed.fastq: fastq that did not have pairs after QC.
    - [`*.stats.txt`](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/samp3.stats.txt): Summary file with numbers of reads before and after QC.
- `ballgown`
  `ballgown` folder. The folder is to be read by `R` package `ballgown` for finding significantly expressed genes. There is one folder per sample.

- `*merged_transcript.gtf`: Non-redundant list of transcripts in GTF format merged from all samples.

- `featureCounts`: A folder containing tables of counts from `featureCounts`.
  - [CDS.count](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/CDS.count):Reads mapped to regions annotated as CDS.
  - [CDS.count.summary](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/CDS.count.summary): Summary of reads mapped and unmapped to CDS.
  - [exon.count](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/exon.count)
  - [exon.count.summary](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/exon.count.summary)
  - [prok_CDS.count](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/prok_CDS.count) : When used `both` option, prokaryote counts are in this file. Eukaryotes are found in file named `euk_CDS.count`
  - [prok_CDS.count.summary](https://raw.githubusercontent.com/mshakya/PyPiReT/master/examples/prok_CDS.count.summary): Corresponding summary file.
- `edgeR`: A folder containing tables and figures processed mainly using R package `edgeR` to detect significantly expressed genes. Based on the options picked, the folder will have either one or two folders, `prokarya` and `eukarya`. Withing these folders there are following files and figures.
  - `*RPKM.csv`: A table with RPKM values for all genes across all samples.
  - `*CPM.csv`: A table with CPM values for all features across all samples
  - `*feature_count_heatmap.pdf`: Heatmap based on count data for the features listed in gff files.
  - `*feature_count_CPM_histogram.pdf`: A histogram of CPMs.
  - `*MDS.pdf`: A MDS plot based on reads mapped to samples.
  - `group1__group2__gene__et.csv`: table with gene name, logFC, logCPM, PValue, and FDR comparing group1 vs. group 2. This one contains all genes that have any counts.
  - `group1__group2__gene__sig.csv`: A subset of `group1__group2__gene__et.csv` with all only genes that are significant based on the specified P-value.


## 4.0 Removing PiReT

For removal, since all dependencies that are not in your system are installed in `PiReT`, delete (`rm -rf`) `PiReT` folder is sufficient to uninstall the package. **Before removing check if your project files are within `PiReT` directory**.


## 5.0 Contributions
- Migun Shakya

## 6.0 Citations:
If you use PiReT please cite following papers:

- **samtools**: Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- **bowtie2**: Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359. [PMID: 22388286]
- **bwa**: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- **DESeq2**: Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. [PMID: 25516281]
- **edgeR**: McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), pp. -9. [PMID: 22287627]
- **HTSeq**: Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq–a Python framework to work with high-throughput sequencing data. Bioinformatics. [PMID: 25260700]
- **hisat2**: Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360. [PMID: 25751142]
- **BEDTools**: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. [PMID: 20110278]
- **GAGE**: Luo, Weijun, Michael S. Friedman, Kerby Shedden, Kurt D. Hankenson, and Peter J. Woolf. 2009. “GAGE: Generally Applicable Gene Set Enrichment for Pathway Analysis.” BMC Bioinformatics 10 (May): 161.
- **Pathview**: Luo, Weijun, and Cory Brouwer. 2013. “Pathview: An R/Bioconductor Package for Pathway-Based Data Integration and Visualization.” Bioinformatics  29 (14). Oxford University Press: 1830–31.
- **Ballgown**: Frazee, Alyssa C., Geo Pertea, Andrew E. Jaffe, Ben Langmead, Steven L. Salzberg, and Jeffrey T. Leek. 2015. “Ballgown Bridges the Gap between Transcriptome Assembly and Expression Analysis.” Nature Biotechnology 33 (3): 243–46.
- **featureCounts**: Liao, Yang, Gordon K. Smyth, and Wei Shi. 2014. “featureCounts: An Efficient General Purpose Program for Assigning Sequence Reads to Genomic Features.” Bioinformatics  30 (7): 923–30.
- **StringTie**: Pertea, Mihaela, Geo M. Pertea, Corina M. Antonescu, Tsung-Cheng Chang, Joshua T. Mendell, and Steven L. Salzberg. 2015. “StringTie Enables Improved Reconstruction of a Transcriptome from RNA-Seq Reads.” Nature Biotechnology 33 (3): 290–95.


## Copyright

Copyright (XXXX).  Triad National Security, LLC. All rights reserved.
 
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National 
Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National 
Nuclear Security Administration.
 
All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National 
Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, 
paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to 
the public, perform publicly and display publicly, and to permit others to do so.

This is open source software; you can redistribute it and/or modify it under the terms of the GPLv3 License. If software 
is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with 
the version available from LANL. Full text of the [GPLv3 License](https://github.com/losalamos/edge/blob/master/LICENSE) can be found in the License file in the main development 
branch of the repository.