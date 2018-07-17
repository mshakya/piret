#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
import glob
from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter, DictParameter, ListParameter, IntParameter
from luigi.util import inherits, requires
from itertools import chain
from plumbum.cmd import FaQCs, cat
import pandas as pd


class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class PairedRunQC(luigi.Task):
    """Running FaQCs."""

    fastqs = ListParameter()
    sample = Parameter()
    num_cpus = IntParameter()
    qc_outdir = Parameter()
    faqc_min_L = IntParameter()
    n_cutoff = IntParameter()

    def requires(self):
        """Require pair of fastq."""
        if isinstance(self.fastqs, (list, tuple)):
            for fqs in self.fastqs:
                fqs_list = fqs.split(",")
                for f_q in fqs_list:
                    return RefFile(f_q)
        elif isinstance(self.fastqs, str):
            return RefFile(self.fastqs.split(":")[0])

    def output(self):
        """QC output."""
        out_file = self.qc_outdir + "/" + self.sample + ".stats.txt"
        return LocalTarget(out_file)

    def run(self):
        """Run the FaQC script."""
        faqc_options = ["-min_L", self.faqc_min_L,
                        "-n", self.n_cutoff,
                        "-t", self.num_cpus,
                        "-prefix", self.sample,
                        "-d", os.path.abspath(self.qc_outdir),
                        "-1", self.fastqs[0],
                        "-2", self.fastqs[1]]
        faqc_cmd = FaQCs[faqc_options]
        faqc_cmd()



class RunAllQC(luigi.WrapperTask):
    """Run all QC."""

    fastq_dic = DictParameter()
    workdir = Parameter()
    num_cpus = IntParameter()
    faqc_min_L = IntParameter()
    n_cutoff = IntParameter()

    def requires(self):
        """A wrapper for running the QC."""
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            if os.path.isdir(trim_dir) is False:
                os.makedirs(trim_dir)
            if isinstance(fastq, (list, tuple)):
                fqs = [fq.replace(';', ',') for fq in fastq]
                i = 1
                for fq in fqs:
                    fq_list = fq.split(",")
                    cp_fq = trim_dir + "/" + samp + "_R" + str(i) + ".fastq"
                    cat_cmd = (cat[fq_list] > cp_fq)
                    cat_cmd()
                    i = i + 1
                yield PairedRunQC(fastqs=[trim_dir + "/" + samp +
                                          "_R1.fastq", trim_dir + "/" +
                                          samp + "_R2.fastq"],
                                  sample=samp,
                                  num_cpus=self.num_cpus,
                                  qc_outdir=trim_dir,
                                  faqc_min_L=self.faqc_min_L,
                                  n_cutoff=self.n_cutoff)

            else:
                if os.path.isdir(trim_dir) is False:
                    os.makedirs(trim_dir)
                fqs = [os.path.abspath(fq) for fq in fastq.split(":")]
                yield PairedRunQC(fastqs=fqs,
                                  sample=samp,
                                  num_cpus=self.num_cpus,
                                  qc_outdir=trim_dir,
                                  faqc_min_L=self.faqc_min_L,
                                  n_cutoff=self.n_cutoff)

@requires(RunAllQC)
class SummarizeQC(luigi.Task):
    """Summarizes FaQC results of all samples into a table"""

    def output(self):
        """QC output."""
        out_file = self.workdir + "/" + "QCsummary.csv"
        return LocalTarget(out_file)

    def run(self):
        """Parse the FaQC stats."""
        summ_dic = {}

        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            filename = trim_dir + "/" + samp + ".stats.txt"
            with open(filename, 'r') as file:
                lines = file.readlines()
                reads_before_trimming = lines[1].split(":")[1].strip()
                read_length = lines[3].split(":")[1].strip()
                reads_aft_trim = lines[9].split(":")[1].split("(")[0].strip()
                summ_dic[samp] = [read_length,
                                  reads_before_trimming,
                                  reads_aft_trim]
        summ_table = pd.DataFrame.from_dict(summ_dic, orient='index')
        summ_table.columns = ["Read Length", "Raw reads", "Reads after QC"]
        out_file = self.workdir + "/" + "QCsummary.csv"
        summ_table.to_csv(out_file)
