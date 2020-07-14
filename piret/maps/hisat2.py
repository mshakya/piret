
"""Luigi Tasks to perform various RNA seq functions

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
"""

from piret.miscs import RefFile
from piret.qc import FaQC
from collections import Counter
from plumbum.cmd import STAR
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import hisat2
from luigi.util import inherits, requires
from luigi import Parameter, IntParameter, DictParameter, ListParameter
from luigi import LocalTarget
from luigi.contrib.external_program import ExternalProgramTask
import os
import sys
import luigi
import pandas as pd
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)


class HisatIndex(ExternalProgramTask):
    """Create Hisat Indices from given fasta file.

    Note: Still using ExternalProgramtask as importing commands with -
    creates error in plumbum

    And, it automatically prints command in log file under INFO
    """

    fasta = Parameter()
    hi_index = Parameter()
    num_cpus = IntParameter()

    def requires(self):
        """Require reference fasta format file."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [RefFile(fa) for fa in fas]
        else:
            return [RefFile(self.fasta)]

    def output(self):
        """Expected index output."""
        hisat_index_ht8l = self.hi_index + ".8.ht2"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["hisat2-build", "-q", "-p", self.num_cpus, self.fasta,
                self.hi_index]


class Hisat(luigi.Task):
    """Mapping the QCed sequences to reference."""

    fastqs = ListParameter()
    indexfile = Parameter()
    workdir = Parameter()
    num_cpus = IntParameter()
    sample = Parameter()
    min_introlen = luigi.IntParameter()
    max_introlen = luigi.IntParameter()
    rna_strandness = Parameter()
    kingdom = luigi.Parameter()

    def requires(self):
        """Check if fastq exist"""
        return [RefFile(os.path.abspath(self.fastqs[0]))]

    def output(self):
        """SAM file output of the mapping."""
        outbam = os.path.join(self.workdir, "processes",
                              "mapping", self.sample, self.sample + ".bam")
        return luigi.LocalTarget(outbam)

    def run(self):
        """Run hisat2."""
        map_dir = os.path.join(self.workdir, "processes",
                               "mapping", self.sample)
        if os.path.isdir(map_dir) is False:
            os.makedirs(map_dir)
        outsam = os.path.join(map_dir, self.sample + ".sam")
        if self.kingdom == "prokarya":
            hisat2_nosplice_option = ["-p", self.num_cpus,
                                      "-x", self.indexfile,
                                      "-1", self.fastqs[0],
                                      "-2", self.fastqs[1],
                                      "-S", outsam,
                                      "--min-intronlen", self.min_introlen,
                                      "--max-intronlen", self.max_introlen,
                                      "--rna-strandness", self.rna_strandness,
                                      "--no-spliced-alignment",
                                      "--no-unal",
                                      "--un-conc",
                                      os.path.join(map_dir,
                                                   "unaligned.fastq"),
                                      "2>", os.path.join(map_dir,
                                                         "mapping.log")]
            hisat2_cmd = hisat2[hisat2_nosplice_option]
            hisat2_cmd()
            self.sam2bam(outsam)
            self.sort_bam(outsam)
        else:
            h2_splice_option = ["-p", self.num_cpus,
                                "-x", self.indexfile,
                                "-1", self.fastqs[0],
                                "-2", self.fastqs[1],
                                "-S", outsam,
                                "--min-intronlen", self.min_introlen,
                                "--max-intronlen", self.max_introlen,
                                "--rna-strandness", self.rna_strandness,
                                "--no-unal",
                                "--un-conc",
                                os.path.join(map_dir,
                                             "unaligned.fastq"),
                                "2>", os.path.join(map_dir,
                                                   "mapping.log")]
            hisat2_cmd = hisat2[h2_splice_option]
            hisat2_cmd()
            self.sam2bam(outsam)
            self.sort_bam(outsam)

    def sam2bam(self, outsam):
        """Convert SAM to BAM file."""
        bam_file = outsam.split(".sam")[0] + ".bam"
        options = ["view", "-bS", "-F",
                   "4", outsam, "-o", bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()

    def sort_bam(self, outsam):
        """Sort BAM file."""
        bam_file = outsam.split(".sam")[0] + ".bam"
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        options = ["sort", bam_file,
                   "-o", sorted_bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()


class HisatMapW(luigi.WrapperTask):
    """A wrapper task for mapping."""

    fastq_dic = DictParameter()
    indexfile = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()
    kingdom = luigi.Parameter()

    def requires(self):
        """A wrapper task for running mapping."""
        splist = [self.workdir + "/" +
                  f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        elif len(splist) == 1:
            splice_file = splist[0]
        else:
            splice_file = ''
        for samp, fastq in self.fastq_dic.items():
            trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
            yield Hisat(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                trim_dir + "/" + samp + ".2.trimmed.fastq"],
                        kingdom=self.kingdom,
                        num_cpus=self.num_cpus,
                        indexfile=self.indexfile,
                        sample=samp,
                        workdir=self.workdir)


@requires(HisatMapW)
class SummarizeHisatMap(luigi.Task):
    """Summarizes mapping results of all samples into a table"""

    def output(self):
        """Mapping Summary Output."""
        out_file = os.path.join(self.workdir, "processes", "mapping",
                                "MapSummary.csv")
        return luigi.LocalTarget(out_file)

    def run(self):
        """Parse the mapping stats."""
        summ_dic = {}
        for samp, fastq in self.fastq_dic.items():
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            filename = map_dir + "/" + "mapping.log"
            with open(filename, 'r') as file:
                lines = file.readlines()
                total_reads = lines[0].split("reads")[0].strip()
                con_unaligned = lines[2].split("(")[0].strip()
                con_aligned = lines[3].split("(")[0].strip()
                multi_aligned = lines[4].split("(")[0].strip()
                summ_dic[samp] = [total_reads,
                                  con_unaligned,
                                  con_aligned,
                                  multi_aligned]
        sm_tbl = pd.DataFrame.from_dict(summ_dic, orient='index')
        sm_tbl.columns = ["pr_re", "concord_unaln",
                          "Concordantly aligned", "Multi aligned"]
        sm_tbl = sm_tbl.astype('int32')
        sm_tbl["perc_unaln"] = (sm_tbl["concord_unaln"]/sm_tbl["pr_re"])*100
        sm_tbl["perc_aln"] = 100-sm_tbl["perc_unaln"]
        out_file = os.path.join(self.workdir, "processes", "mapping",
                                "MapSummary.csv")
        sm_tbl = sm_tbl.round(2)
        sm_tbl = sm_tbl.rename(columns={'pr_re': 'Paired Reads',
                                        'concord_unaln': 'Concordantly Unaligned',
                                        'perc_unaln': "Percentage Unaligned",
                                        'perc_aln': "Percentage Aligned"})
        sm_tbl.to_csv(out_file)
