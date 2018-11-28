#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts and stringtie
"""

from __future__ import print_function
import os
import luigi
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from luigi.contrib.external_program import ExternalProgramTask
from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter, IntParameter, DictParameter, ListParameter
from luigi.util import inherits, requires
import subprocess
from pypiret import FaQC
from plumbum.cmd import gffread, hisat2
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import STAR
import pandas as pd
from sys import stderr, exit
import logging
from collections import defaultdict as dd, Counter
# from Bio import SeqIO

class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class GetChromName(luigi.Task):
    """Extract chromosomes of euk and prok fasta."""

    prok_ref = Parameter()
    euk_ref = Parameter()
    workdir = Parameter()

    def requires(self):
        """Require slist of rerence sequence to be made."""
        return [RefFile(self.euk_ref), RefFile(self.prok_ref)]

    def output(self):
        """Split Files."""
        prok_chrom_file = self.workdir + "/" + "prok.chroms"
        euk_chrom_file = self.workdir + "/" + "euk.chroms"
        return [luigi.LocalTarget(prok_chrom_file),
                luigi.LocalTarget(euk_chrom_file)]

    def run(self):
        """Run the subprocess for grep and sed."""
        prok_chrom_file = self.workdir + "/" + "prok.chroms"
        prok_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.prok_ref,
                                                                          prok_chrom_file)
        subprocess.Popen(prok_cmd, shell=True)

        euk_chrom_file = self.workdir + "/" + "euk.chroms"
        euk_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.euk_ref,
                                                                         euk_chrom_file)
        subprocess.Popen(euk_cmd, shell=True)

class HisatIndex(ExternalProgramTask):
    """Create Hisat Indices from given fasta file."""

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
        hisat_index_ht8l = self.hi_index + ".8.ht2l"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["hisat2-build",
                "--large-index",
                "-q", "-p", self.num_cpus,
                self.fasta, self.hi_index]


class SAMindex(luigi.Task):
    """Create Hisat Indices from given fasta file."""

    fasta = Parameter()
    workdir = Parameter()

    def output(self):
        """Require reference fasta format file."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [LocalTarget(self.workdir + "/" + os.path.basename(fa) +
                                ".bedfile") for fa in fas]
        else:
            return [LocalTarget(self.workdir + "/" +
                                os.path.basename(self.fasta) +
                                ".bedfile")]

    def make_index(self, ref, workdir):
        """A function to make index from sam files."""
        index_options = ["faidx", ref]
        mv_options = [ref + ".fai", workdir]
        samtools_cmd = samtools[index_options]
        mv_cmd = mv[mv_options]
        samtools_cmd()
        mv_cmd()
        fa_name = os.path.basename(ref)
        return os.path.abspath(self.workdir + "/" + fa_name + ".fai")

    def create_bedfile(self, index_file):
        """Makes bed file from the sam index file."""
        out_bedfile = index_file.split(".fai")[0] + ".bedfile"
        bedfile_opt = ['BEGIN {FS="\t"}; {print $1 FS "0" FS $2}',
                       index_file]
        awk_cmd = ((awk[bedfile_opt]) > out_bedfile)
        awk_cmd()

    def requires(self):
        """Expected index output."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [RefFile(fa) for fa in fas]
        else:
            return [RefFile(self.fasta)]

    def run(self):
        """Run hisat2-build command."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            for fa in fas:
                ind_file = self.make_index(fa, self.workdir)
                self.create_bedfile(ind_file)


class GFF2GTF(luigi.Task):
    """Converts GFF to GTF format."""

    gff_file = Parameter()
    workdir = Parameter()

    def requires(self):
        """Require reference gff(s) file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [RefFile(os.path.abspath(gff)) for gff in gffs]
        else:
            return [RefFile(os.path.abspath(self.gff_file))]

    def output(self):
        """Output converted gtf file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                os.path.abspath(gff).split("/")[-1].rsplit(".", 1)[0] +
                                ".gtf") for gff in gffs]
        else:
            gff = os.path.abspath(self.gff_file)
            return LocalTarget(self.workdir + "/" +
                               gff.split("/")[-1].rsplit(".", 1)[0] +
                               ".gtf")

    def run(self):
        """Main."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                gff_abs = os.path.abspath(gff)
                out_file = self.workdir + "/" +\
                    gff_abs.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
                gffread_option = [gff_abs, "-O", "-o", out_file]
                gffread_cmd = gffread[gffread_option]
                gffread_cmd()
        else:
            gff = os.path.abspath(self.gff_file)
            out_file = self.workdir + "/" +\
                gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
            gffread_option = [self.gff_file, "-O","-o", out_file]
            gffread_cmd = gffread[gffread_option]
            gffread_cmd()


@requires(GFF2GTF)
class CreateSplice(ExternalProgramTask):
    """Find splice sites off gtf file."""

    def output(self):
        """Splice site output."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                os.path.abspath(gff).split("/")[-1].split(".")[0] +
                                ".splice") for gff in gffs]
        else:
            gff = os.path.abspath(self.gff_file)
            return LocalTarget(self.workdir + "/" +
                               gff.split("/")[-1].split(".")[0] +
                               ".splice")

    def extract_splice_sites(self, gtf_file, splice, verbose=False):
        """"Function from hisat2_extract_splice_sites.py script."""
        genes = dd(list)
        trans = {}

        # Parse valid exon lines from the GTF file into a dict by transcript_id
        for line in gtf_file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                chrom, source, feature, left, right, score, \
                    strand, frame, values = line.split('\t')
            except ValueError:
                continue
            left, right = int(left), int(right)

            if feature != 'exon' or left >= right:
                continue

            values_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

            if 'gene_id' not in values_dict or \
                    'transcript_id' not in values_dict:
                continue

            transcript_id = values_dict['transcript_id']
            if transcript_id not in trans:
                trans[transcript_id] = [chrom, strand, [[left, right]]]
                genes[values_dict['gene_id']].append(transcript_id)
            else:
                trans[transcript_id][2].append([left, right])

        # Sort exons and merge where separating introns are <=5 bps
        for tran, [chrom, strand, exons] in trans.items():
                exons.sort()
                tmp_exons = [exons[0]]
                for i in range(1, len(exons)):
                    if exons[i][0] - tmp_exons[-1][1] <= 5:
                        tmp_exons[-1][1] = exons[i][1]
                    else:
                        tmp_exons.append(exons[i])
                trans[tran] = [chrom, strand, tmp_exons]

        # Calculate and print the unique junctions
        junctions = set()
        for chrom, strand, exons in trans.values():
            for i in range(1, len(exons)):
                junctions.add((chrom, exons[i-1][1], exons[i][0], strand))
        junctions = sorted(junctions)
        with open(splice, 'w') as out:
            for chrom, left, right, strand in junctions:
                # Zero-based offset
                out.write('{}\t{}\t{}\t{}\n'.format(chrom, left-1, right-1, strand))
                # print('{}\t{}\t{}\t{}'.format(chrom, left-1, right-1, strand))

        # Print some stats if asked
        if verbose:
            exon_lengths, intron_lengths, trans_lengths = \
                Counter(), Counter(), Counter()
            for chrom, strand, exons in trans.values():
                tran_len = 0
                for i, exon in enumerate(exons):
                    exon_len = exon[1]-exon[0]+1
                    exon_lengths[exon_len] += 1
                    tran_len += exon_len
                    if i == 0:
                        continue
                    intron_lengths[exon[0] - exons[i-1][1]] += 1
                trans_lengths[tran_len] += 1

            print('genes: {}, genes with multiple isoforms: {}'.format(
                  len(genes), sum(len(v) > 1 for v in genes.values())),
                  file=stderr)
            print('transcripts: {}, transcript avg. length: {:d}'.format(
                    len(trans), sum(trans_lengths.elements())/len(trans)),
                  file=stderr)
            print('exons: {}, exon avg. length: {:d}'.format(
                    sum(exon_lengths.values()),
                    sum(exon_lengths.elements())/sum(exon_lengths.values())),
                  file=stderr)
            print('introns: {}, intron avg. length: {:d}'.format(
                    sum(intron_lengths.values()),
                    sum(intron_lengths.elements())/sum(intron_lengths.values())),
                  file=stderr)
            print('average number of exons per transcript: {:d}'.format(
                    sum(exon_lengths.values())/len(trans)),
                  file=stderr)


    def run(self):
        """Main."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                gtf_file = self.workdir + "/" +\
                    gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
                out_file = self.workdir + "/" +\
                    gff.split("/")[-1].rsplit(".", 1)[0] + ".splice"
                self.extract_splice_sites(gtf_file, out_file)
        else:
            gff = os.path.abspath(self.gff_file)
            gtf_file = self.workdir + "/" +\
                gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
            out_file = self.workdir + "/" +\
                self.gff_file.split("/")[-1].rsplit(".", 1)[0] + ".splice"
            self.extract_splice_sites(gtf_file, out_file)


# @requires(FaQC.PairedRunQC)
class Hisat(luigi.Task):
    """Mapping the QCed sequences to reference."""

    fastqs = ListParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    outsam = Parameter()
    map_dir = Parameter()
    workdir = Parameter()
    num_cpus = Parameter()
    sample=Parameter()

    def output(self):
        """SAM file output of the mapping."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        return luigi.LocalTarget(bam_file)

    def run(self):
        """Run hisat2."""
        if self.spliceFile is "":
            hisat2_nosplice_option = ["-p", self.num_cpus,
                                      "-x", self.indexfile,
                                      "-1", self.fastqs[0],
                                      "-2", self.fastqs[1],
                                      "-S", self.outsam,
                                      "--no-unal",
                                      "--un-conc",
                                      os.path.join(self.map_dir,
                                                   "unaligned.fastq"),
                                      "2>", os.path.join(self.map_dir,
                                                         "mapping.log")]
            hisat2_cmd = hisat2[hisat2_nosplice_option]
            hisat2_cmd()
            self.sam2bam()
            self.sort_bam()
        else:
            h2_splice_option = ["--known-splicesite-infile", self.spliceFile,
                                "-p", self.num_cpus,
                                "-x", self.indexfile,
                                "-1", self.fastqs[0],
                                "-2", self.fastqs[1],
                                "-S", self.outsam,
                                "--no-unal",
                                "--un-conc",
                                os.path.join(self.map_dir,
                                             "unaligned.fastq"),
                                "2>", os.path.join(self.map_dir,
                                                   "mapping.log")]
            hisat2_cmd = hisat2[h2_splice_option]
            hisat2_cmd()
            self.sam2bam()
            self.sort_bam()

    def sam2bam(self):
        """Convert SAM to BAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        options = ["view", "-bS", "-F",
                   "4", self.outsam, "-o", bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()

    def sort_bam(self):
        """Sort BAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
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
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield Hisat(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                trim_dir + "/" + samp + ".2.trimmed.fastq"],
                        num_cpus=self.num_cpus,
                        indexfile=self.indexfile,
                        spliceFile=splice_file,
                        outsam=map_dir + "/" + samp + ".sam",
                        map_dir=map_dir,
                        sample=samp,
                        workdir=self.workdir)


@requires(HisatMapW)
class SummarizeHisatMap(luigi.Task):
    """Summarizes mapping results of all samples into a table"""

    def output(self):
        """Mapping Summary Output."""
        out_file = self.workdir + "/" + "MapSummary.csv"
        return luigi.LocalTarget(out_file)

    def run(self):
        """Parse the FaQC stats."""
        summ_dic = {}
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
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
        summ_table = pd.DataFrame.from_dict(summ_dic, orient='index')
        summ_table.columns = ["Paired reads", "Concordantly unaligned",
                              "Concordantly aligned", "Multi aligned"]
        out_file = self.workdir + "/" + "MapSummary.csv"
        summ_table.to_csv(out_file)

class STARindex(luigi.Task):
    """ Creates STAR index."""
    fasta = Parameter()
    num_cpus = IntParameter()
    gff_file = Parameter()
    stardb_dir= Parameter()

    def requires(self):
        return[RefFile(self.fasta)]
    
    def output(self):
        """Expected index output"""
        return LocalTarget(os.path.join(self.stardb_dir, "chrName.txt"))
    
    def run(self):
        if os.path.exists(self.stardb_dir) is False:
            os.makedirs(self.stardb_dir)
        ind_opt = ["--runMode", "genomeGenerate",
                   "--runThreadN", self.num_cpus,
                   "--genomeDir", self.stardb_dir,
                   "--genomeFastaFiles", self.fasta,
                   "--sjdbGTFfile", self.gff_file]
        star_cmd = STAR[ind_opt]
        logger = logging.getLogger('luigi-interface')
        logger.info(star_cmd)
        star_cmd()


class map_star(luigi.Task):
    """Mapping the QCed sequences to reference using star aligner."""
    fastqs = ListParameter()
    stardb_dir = Parameter()
    map_dir = Parameter()
    sample = Parameter()
    num_cpus = IntParameter()

    # def requires(self):
    #     """See if input file exist."""
    #     return [luigi.LocalTarget(self.fastqs[0]),
    #             luigi.LocalTarget(self.fastqs[1])]

    def output(self):
        """SAM file output of the mapping."""
        bam_file = os.path.join(self.map_dir, self.sample) + "_srt.bam"
        return luigi.LocalTarget(bam_file)

    def run(self):
        """Run star"""
        if os.path.exists(self.map_dir) is False:
            os.makedirs(self.map_dir)
        star_option = ["--genomeDir", self.stardb_dir,
                       "--runThreadN", self.num_cpus,
                       "--outFileNamePrefix", os.path.join(self.map_dir, self.sample + "_"),
                       "--outSAMtype", "BAM", "SortedByCoordinate",
                       "--readFilesIn", self.fastqs[0], self.fastqs[1]]
        star_cmd = STAR[star_option]
        star_cmd()
        bam_file = os.path.join(self.map_dir, self.sample) + "_Aligned.sortedByCoord.out.bam"
        mv_list = [bam_file, os.path.join(self.map_dir, self.sample) + "_srt.bam"]
        mv_cmd = mv[mv_list]
        mv_cmd()


class map_starW(luigi.WrapperTask):
    """A wrapper task for mapping using star."""

    fastq_dic = luigi.DictParameter()
    stardb_dir = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()

    def requires(self):
        """A wrapper task for mapping using STAR."""
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield map_star(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                   trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                   stardb_dir=self.stardb_dir, map_dir=map_dir,
                                   sample=samp, num_cpus=self.num_cpus)

@requires(map_starW)
class SummarizeStarMap(luigi.Task):
    """Summarizes mapping results of all samples into a table"""

    def output(self):
        """Mapping Summary Output."""
        out_file = self.workdir + "/" + "MapSummary.csv"
        return luigi.LocalTarget(out_file)

    def run(self):
        """Parse the mapping stats."""
        summ_dic = {}
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            filename = map_dir + "/" + samp + "_Log.final.out"
            with open(filename, 'r') as file:
                lines = file.readlines()
                total_reads = lines[5].split("|")[1].strip()
                unq_map_reads = lines[8].split("|")[1].strip()
                multi_map_reads = lines[23].split("|")[1].strip()
                multiX_map_reads = lines[25].split("|")[1].strip()
                summ_dic[samp] = [total_reads,
                                  unq_map_reads,
                                  multi_map_reads,
                                  multiX_map_reads]
        summ_table = pd.DataFrame.from_dict(summ_dic, orient='index')
        summ_table.columns = ["Paired reads", "Uniquely mapped reads",
                              "mapped to multiple loci", "mapped to too many loci"]
        out_file = self.workdir + "/" + "MapSummary.csv"
        summ_table.to_csv(out_file)

@requires(Hisat)
class RefNames(luigi.Task):
    """Extract the name of chromosomes where reads were mapped in BAM file."""

    map_ref = Parameter()

    def output(self):
        """Split Files."""
        return luigi.LocalTarget(self.map_ref)

    def run(self):
        """Extract names of chromosome from SAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        grep_cmd = "samtools view -h %s |\
                    grep -v -e '^.SQ' |\
                    grep -v -e 'PN:' |\
                    grep -v -e 'SO:' |\
                    awk -F '\t' '{print $3}' |\
                    sort | uniq |\
                    grep -v '*' > %s" % (sorted_bam_file, self.map_ref)
        subprocess.Popen(grep_cmd, shell=True)


@inherits(HisatMapW)
class GetRefNames(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splist = [self.workdir + "/" +
                  f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        else:
            splice_file = splist[0]
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield RefNames(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                   trim_dir + "/" + samp + ".2.trimmed.fastq"],
                           num_cpus=self.num_cpus,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           outsam=map_dir + "/" + samp + ".sam",
                           map_ref=map_dir + "/" + samp + ".reflist",
                           ref_file=self.ref_file,
                           sample=samp,
                           qc_outdir=trim_dir,
                           map_dir=map_dir)


@inherits(GFF2GTF)
@requires(Hisat)
class StringTieScores(luigi.Task):
    """Calculate scores using string tie."""

    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()

    def output(self):
        """Index output."""
        return LocalTarget(self.out_gtf)

    def run(self):
        """Run stringtie."""
        stringtie_opt = ["-o", self.out_gtf,
                         "-p", self.num_cpus,
                         "-G", self.gff_file,
                         "-C", self.out_cover,
                         "-A", self.out_abun,
                         self.in_bam_file]
        stringtie_cmd = stringtie[stringtie_opt]
        stringtie_cmd()


@inherits(HisatMapW)
class StringTieScoresW(luigi.WrapperTask):
    """Wrapper function for stringtie in all samples"""

    gff_file = Parameter()
    kingdom = Parameter()

    def requires(self):
        """A wrapper for running Stringtie scores on all samples."""
        splst = [self.workdir + "/" +
                 f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splst) > 1:
            splice_file = ','.join(splst)
        elif len(splst) == 1:
            splice_file = splst[0]
        else:
            splice_file = ""
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            stng_dir = self.workdir + "/" + samp + "/stie_results"
            if os.path.isdir(stng_dir) is False:
                os.makedirs(stng_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                if self.kingdom == 'prokarya':
                    apd = '_prok'
                elif self.kingdom == 'eukarya':
                    apd = '_euk'

                yield StringTieScores(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                              trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      gff_file=self.gff_file,
                                      out_gtf=stng_dir + "/" + samp +  apd + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp +  apd + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp +  apd + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                      workdir=self.workdir,
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)
            elif self.kingdom == 'both':
                prok_gff = self.gff_file.split(",")[0]
                euk_gff = self.gff_file.split(",")[1]
                yield StringTieScores(fastqs=[trim_dir + "/" + samp +
                                              ".1.trimmed.fastq",
                                              trim_dir + "/" + samp +
                                              ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      out_gtf=stng_dir + "/" + samp + "_prok" + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp + "_prok" + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp + "_prok" + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt_prok.bam",
                                      workdir=self.workdir,
                                      gff_file=prok_gff,
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)
                yield StringTieScores(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                      trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      out_gtf=os.path.join(stng_dir, samp + "_euk_sTie.gtf"),
                                      out_cover=os.path.join(stng_dir, samp + "_euk_covered_sTie.gtf"),
                                      out_abun=os.path.join(stng_dir, samp + "_euk_sTie.tab"),
                                      in_bam_file=map_dir + "/" + samp + "_srt_prok.bam",
                                      workdir=self.workdir,
                                      gff_file=euk_gff,
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)

@requires(Hisat)
class Split2ProkEuk(luigi.Task):
    """ Split BAM file to prok and euk"""

    ref_file = Parameter()
    workdir = Parameter()

    def output(self):
        """Split Files."""
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        return [luigi.LocalTarget(prok_out), luigi.LocalTarget(euk_out)]

    def split_aln_file(self, workdir, fasta, bam_file, out_bamfile):
        """Split"""
        bed_file = workdir + "/" + os.path.basename(fasta) + ".bedfile"
        samtools_opt = ["view", "-b", "-L", bed_file, bam_file]
        samtools_cmd = (samtools[samtools_opt]) > out_bamfile
        samtools_cmd()

    def run(self):
        """ Split BAM file to proks and euks"""
        fastas = self.ref_file.split(",")
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        self.split_aln_file(self.workdir, fastas[0], bam_file, prok_out)
        self.split_aln_file(self.workdir, fastas[1], bam_file, euk_out)


@inherits(HisatMapW)
class Split2ProkEukW(luigi.WrapperTask):
    """Group chromosomes to prokaryotic and eukaryotic."""
    ref_file = Parameter()
    workdir = Parameter()

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splst = [self.workdir + "/" +
                 f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splst) > 1:
            splice_file = ','.join(splst)
        else:
            splice_file = splst[0]
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield Split2ProkEuk(fastqs=[trim_dir + "/" + samp +
                                        ".1.trimmed.fastq",
                                        trim_dir + "/" + samp +
                                        ".2.trimmed.fastq"],
                                num_cpus=self.num_cpus,
                                indexfile=self.indexfile,
                                spliceFile=splice_file,
                                outsam=map_dir + "/" + samp + ".sam",
                                map_dir=map_dir,
                                ref_file=self.ref_file,
                                
                                workdir=self.workdir,
                                sample=samp,
                                qc_outdir=trim_dir)

