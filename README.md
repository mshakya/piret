<!-- [![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://bioconda.github.io/recipes/viral-ngs/README.html) -->

[![Build Status](https://travis-ci.org/mshakya/PyPiReT.svg?branch=master)](https://travis-ci.org/mshakya/PyPiReT)
[![codecov](https://codecov.io/gh/mshakya/PyPiReT/branch/master/graph/badge.svg)](https://codecov.io/gh/mshakya/PyPiReT)

#PiReT

Pipeline for Reference based Transcriptomics.

#Overview

![alt tag](image/Overview_pipeline.jpg)

## Installing PiReT

Please download PiReT directly from the [github](https://github.com/mshakya/PyPiReT).

```
git clone https://github.com/mshakya/PyPiReT.git
```

`cd` into the `PyPiReT` directory

```
cd PyPiReT
./INSTALL.SH
```

PiReT uses bioinformatic tools, many of which are available in [bioconda](https://bioconda.github.io). For installing `PiReT` we have provided a script `bioconda_INSTALL.sh` that checks for required dependencies (including their versions) are installed and in your path, and installs it in directories within `PiReT` if not found. Additionally, `sudo` privileges are not needed for installation. A log of all installation can be found in `install.log`

##Test
We have provided test data set to check if the installation was successful or not. `fastq` files can be found in `tests/fastqs` and corresponding reference fasta files are found in `tests/data`. To run the test, from within `PyPiReT` directory:

```
# if you are in a LINUX system:
sh ./test_pipeline_linux.sh

# if you are in Mac OS X:
sh ./test_pipeline_MacOSX.sh
```
These shell script automatically creates `experimental_design.txt` and runs the pipeline.


##Dependencies
PiReT requires following dependencies, all of which should be installed and in the PATH. All of the dependencies will be installed by `INSTALL.sh`.

### Programming/Scripting languages
- [Python >=v2.7](https://www.python.org/downloads/release/python-2712/)
    - The pipeline is not compatible with Python v3.0 or higher.
- [Perl >=v5.16.3](https://www.perl.org/get.html)
    - The pipeline has only been tested in v5.16.3 and v5.22.0
- [R >=v3.3.1](https://www.r-project.org)

### Installing dependencies
This is the core list of dependencies. However, there are secondary dependencies for many of the listed tools, which will also be installed by `bioconda`.
- [conda v4.2.13](http://conda.pydata.org/docs/index.html)
    If conda is not installed, `bioconda_INSTALL.sh` will download and install [miniconda](http://conda.pydata.org/miniconda.html), a "mini" version of `conda` that only installs handful of packages compared to [anaconda](https://docs.continuum.io/anaconda/pkg-docs)

### Third party softwares/packages
- [jellyfish (v2.2.6)](http://www.genome.umd.edu/jellyfish.html)
- [samtools (v1.3.1)](http://www.htslib.org)
- [HiSat2 (v2.0.5)](https://ccb.jhu.edu/software/hisat/index.shtml)
- [bedtools (v2.26.0)](http://bedtools.readthedocs.io/en/latest/index.html)
- [gffread (v0.9.6)](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_dl)

### R packages
- [edgeR (v3.14.0)](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [DEseq2 (v1.12.4)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### Python packages
- [luigi](https://github.com/spotify/luigi)
- [numpy (v1.1.12)](http://www.numpy.org)
- [matplotlib (v1.5.3)](http://matplotlib.org)

### Perl modules
- [Parallel::ForkManager (v1.17)](http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm)
- [String::Approx (v3.27)](http://search.cpan.org/dist/String-Approx/Approx.pm)

## Running PiReT


```
    bin/runPiReT -test_kingdom both \
    -significant_pvalue 0.001 -exp experimental_design.txt \
    -d pipeline_test_both \
    -prokaryote_fasta data/test_prok.fa \
    -eukarya_fasta data/eukarya_test.fa -index_ref_bt2 test_index \
    -gff_eukarya data/eukarya_test.gff3 -gff_prokaryote data/test_prok.gff \
    -test_method both -gene_coverage_fasta data/test_prok.fa
```

`-d`: working directory where all output files/directories will be written, users must have write permission.

`-prokaryote_fasta`: comma-separated list of reference genome (prokaryote) fasta files. [optional]

`-gff_prokaryote`: comma-separated list of gff files for corresponding reference genome fasta files (contig names must match reference sequence header). [optional]

`-eukarya_fasta` : comma-separated list of reference genome (eukarya) fasta files. [optional]

`-gff_eukarya`: comma-separated list of gff files for corresponding reference genome fasta files (contig names must match reference sequence header). [optional]

`-index_ref_bt2`: HISAT2 mapping index file, if the file exists, pipeline skips this step. [optional]

`-gene_coverage_fasta`: fasta file  (for directional coverage analysis, sequence  must be part of prokaryote mapping reference sequence). [optional]

`-test_kingdom`: desired differential gene expression analysis (`both` (for both eukarya and prokaryote), `prokaryote`, or `eukarya` (default:`prokaryote`));

`-test_method`: method for determining differentially expressed genes. Options are `EdgeR`, `DeSeq2` (For Deseq2, must have have at least 3 replicates for each group), and `both`. `default`: `both`. 

<!-- `-cpu`: number of CPU to be used (default 1) -->

`-BAM_ready`: if mapping file are provided for samples by users (`yes` or `no`). default: `no`

`-significant_pvalue`: floating number cutoff to define significant differentially express genes, (default=0.001)

`-exp`: A tab delimited file that contains at least 3 columns with following header `ID`, `Rawread_files`, and  `group`. `Rawread_files` must have an absolute path.

`-pair_comparison`: tab delimited file describing pairwise comparison. If the file is not specified, all possible pairwise analysis will be conducted.


## Whats in the working directory (-d)?

Here are the list of directories that will be in `working directory`.

```

ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'

|-differential_gene
   |---prokaryote
   |-----test_prok
   |-------Deseq
   |-------EdgeR
   |-------figures
   |-------significant_gene
   |---------pathway
   |-logdir
   |-samp2
   |---mapping_results
   |---trimming_results
   |-samp3
   |---mapping_results
   |---trimming_results
   |-sum_gene_count
   |---read_count
   |-----prokaryote
   |-------test_prok
   |---tmp_count
   |-----prokaryote
   |-------test_prok

```

`differential_gene`: contains sub-folders with `EdgeR` and `DeSeq` results (when provided) for `prokaryote` and `eukaryote` or `both`. Direct sub-directory of `differential_gene` can be either `prokaryote`, `eukarya`, or both of those. The folder within it are then named based on the given `gff` files corresponding usually to one organism. The file and directory structure within each folder are mostly similar with few differences, all of which are listed and described below.

- `eukarya/splice_sites_gff.txt`: contains known splice sites, generated using `scripts/extract_splice_sites.py`, a python script part of *HISAT*.
- `sum_exp_stats.txt`: Summary table of number of reads after each major processing (filtering and mapping) of files.
- `RPKM_all_gene.txt`: A table of RPKM calculated per features for each samples.
- `reads.table.txt` : A table of reads mapped to features for each samples.
- `readcounts.experiment.txt`: table similar to experimental design file with location of `htseq-count` results full path.

`process.log`: report of all the commands/scripts/ that were ran as part of the pipeline.

`error.log`: any error are reported here.

`samp2`: The name of this directory corresponds to sample name. Within this folder there are two sub-folders:

- `mapping_results`
    This folder contains reads mapped using *HISAT2* in following formats. If `splice_sites_gff.txt` is present, **HISAT2** aligns based on known splice sites (`splice_sites_gff.txt`).
        - `.sam`: outputs of *HISAT2* (`forward`, `backward`, `paired`, `Notproperpaired`) and sorted `.sam` files
        - `.bam`: generated from `.sam` with `samtools view -bt index_file .sam < .bam`
        - `.bedgraph`: bedgraph summaries of feature coverage produced using`genomeCoverageBed -split -bg -ibam`
        - `.mapping.log`: Alignment summary file from `HISAT2`

- 2. `trimming_results`
            This folder contains results of quality trimming or filtering. This folder was generated using the same script that filteres reads in [EDGE](https://bioedge.lanl.gov/edge_ui/) pipeline.


`sum_gene_count`: directory with results of reads count per sample. Calculated using `htseq-count  -t gene -q -i locus_tag`. Also see `readcounts.experiment.txt`.

`eukarya.fai`: Indexed reference sequence from `eukarya.fa` using `samtools faidx`. A four column table with NAME, LENGTH, OFFSET, LINEBASES, and LINEWIDTH 

`process_current.log`: Created at the end of the pipeline indicating, successful run.

## Removing PiReT

For removal, since all dependencies that are not in your system are installed in `PiReT`, delete (`rm -rf`) `PiReT` folder is sufficient to uninstall the package. **Before removing check if your project files are within `PiReT` directory**.


##Contributions
- Migun Shakya
- Shihai Feng

## Citations:
If you use PiReT please cite following papers:

- **samtools**: Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]
- **bowtie2**: Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357-359. [PMID: 22388286]
- **bwa**: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- **DESeq2**: Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. [PMID: 25516281]
- **EdgeR**: McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), pp. -9. [PMID: 22287627]
- **HTSeq**: Anders, S., Pyl, P. T., & Huber, W. (2014). HTSeq–a Python framework to work with high-throughput sequencing data. Bioinformatics. [PMID: 25260700]
- **HISAT2**: Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360. [PMID: 25751142]
- **BEDTools**: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842. [PMID: 20110278]
