#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
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
from plumbum.cmd import hisat2
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import STAR
import pandas as pd
from sys import stderr, exit
import logging
from collections import defaultdict as dd, Counter

class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


# class GetChromName(luigi.Task):
#     """Extract chromosomes of euk and prok fasta."""

#     prok_ref = Parameter()
#     euk_ref = Parameter()
#     workdir = Parameter()

#     def requires(self):
#         """Require slist of rerence sequence to be made."""
#         return [RefFile(self.euk_ref), RefFile(self.prok_ref)]

#     def output(self):
#         """Split Files."""
#         prok_chrom_file = self.workdir + "/" + "prok.chroms"
#         euk_chrom_file = self.workdir + "/" + "euk.chroms"
#         return [luigi.LocalTarget(prok_chrom_file),
#                 luigi.LocalTarget(euk_chrom_file)]

#     def run(self):
#         """Run the subprocess for grep and sed."""
#         prok_chrom_file = self.workdir + "/" + "prok.chroms"
#         prok_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.prok_ref,
#                                                                           prok_chrom_file)
#         subprocess.Popen(prok_cmd, shell=True)

#         euk_chrom_file = self.workdir + "/" + "euk.chroms"
#         euk_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.euk_ref,
#                                                                          euk_chrom_file)
#         subprocess.Popen(euk_cmd, shell=True)


class SAMindex(luigi.Task):
    """Create Hisat Indices from given fasta file."""

    fasta = Parameter()
    workdir = Parameter()

    def output(self):
        """Require reference fasta format file."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [LocalTarget(os.path.join(self.workdir, "processes",
                                             "novel",
                                             os.path.basename(fa) +
                                             ".bedfile")) for fa in fas]

    def make_index(self, ref):
        """A function to make index from sam files."""
        index_options = ["faidx", ref]
        novel_folder = os.path.join(self.workdir, "processes", "novel")
        if os.path.exists(novel_folder) is False:
            os.makedirs(novel_folder)
        mv_options=[ref + ".fai", os.path.join(novel_folder)]
        samtools_cmd=samtools[index_options]
        mv_cmd=mv[mv_options]
        samtools_cmd()
        mv_cmd()
        fa_name=os.path.basename(ref)
        return os.path.join(self.workdir, "processes", "novel",
                            fa_name + ".fai")

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
                ind_file = self.make_index(fa)
                self.create_bedfile(ind_file)


# class GFF2GTF(luigi.Task):
#     """Converts GFF to GTF format."""

#     gff_file = Parameter()
#     workdir = Parameter()

#     def requires(self):
#         """Require reference gff(s) file."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             return [RefFile(os.path.abspath(gff)) for gff in gffs]
#         else:
#             return [RefFile(os.path.abspath(self.gff_file))]

#     def output(self):
#         """Output converted gtf file."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             return [LocalTarget(self.workdir + "/" +
#                                 os.path.abspath(gff).split("/")[-1].rsplit(".", 1)[0] +
#                                 ".gtf") for gff in gffs]
#         else:
#             gff = os.path.abspath(self.gff_file)
#             return LocalTarget(self.workdir + "/" +
#                                gff.split("/")[-1].rsplit(".", 1)[0] +
#                                ".gtf")

#     def run(self):
#         """Main."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             for gff in gffs:
#                 gff_abs = os.path.abspath(gff)
#                 out_file = self.workdir + "/" +\
#                     gff_abs.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
#                 gffread_option = [gff_abs, "-O", "-o", out_file]
#                 gffread_cmd = gffread[gffread_option]
#                 gffread_cmd()
#         else:
#             gff = os.path.abspath(self.gff_file)
#             out_file = self.workdir + "/" +\
#                 gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
#             gffread_option = [self.gff_file, "-O","-o", out_file]
#             gffread_cmd = gffread[gffread_option]
#             gffread_cmd()


# @requires(GFF2GTF)
# class CreateSplice(ExternalProgramTask):
#     """Find splice sites off gtf file."""

#     def output(self):
#         """Splice site output."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             return [LocalTarget(self.workdir + "/" +
#                                 os.path.abspath(gff).split("/")[-1].split(".")[0] +
#                                 ".splice") for gff in gffs]
#         else:
#             gff = os.path.abspath(self.gff_file)
#             return LocalTarget(self.workdir + "/" +
#                                gff.split("/")[-1].split(".")[0] +
#                                ".splice")

#     def extract_splice_sites(self, gtf_file, splice, verbose=False):
#         """"Function from hisat2_extract_splice_sites.py script."""
#         genes = dd(list)
#         trans = {}

#         # Parse valid exon lines from the GTF file into a dict by transcript_id
#         for line in gtf_file:
#             line = line.strip()
#             if not line or line.startswith('#'):
#                 continue
#             if '#' in line:
#                 line = line.split('#')[0].strip()

#             try:
#                 chrom, source, feature, left, right, score, \
#                     strand, frame, values = line.split('\t')
#             except ValueError:
#                 continue
#             left, right = int(left), int(right)

#             if feature != 'exon' or left >= right:
#                 continue

#             values_dict = {}
#             for attr in values.split(';')[:-1]:
#                 attr, _, val = attr.strip().partition(' ')
#                 values_dict[attr] = val.strip('"')

#             if 'gene_id' not in values_dict or \
#                     'transcript_id' not in values_dict:
#                 continue

#             transcript_id = values_dict['transcript_id']
#             if transcript_id not in trans:
#                 trans[transcript_id] = [chrom, strand, [[left, right]]]
#                 genes[values_dict['gene_id']].append(transcript_id)
#             else:
#                 trans[transcript_id][2].append([left, right])

#         # Sort exons and merge where separating introns are <=5 bps
#         for tran, [chrom, strand, exons] in trans.items():
#                 exons.sort()
#                 tmp_exons = [exons[0]]
#                 for i in range(1, len(exons)):
#                     if exons[i][0] - tmp_exons[-1][1] <= 5:
#                         tmp_exons[-1][1] = exons[i][1]
#                     else:
#                         tmp_exons.append(exons[i])
#                 trans[tran] = [chrom, strand, tmp_exons]

#         # Calculate and print the unique junctions
#         junctions = set()
#         for chrom, strand, exons in trans.values():
#             for i in range(1, len(exons)):
#                 junctions.add((chrom, exons[i-1][1], exons[i][0], strand))
#         junctions = sorted(junctions)
#         with open(splice, 'w') as out:
#             for chrom, left, right, strand in junctions:
#                 # Zero-based offset
#                 out.write('{}\t{}\t{}\t{}\n'.format(chrom, left-1, right-1, strand))
#                 # print('{}\t{}\t{}\t{}'.format(chrom, left-1, right-1, strand))

#         # Print some stats if asked
#         if verbose:
#             exon_lengths, intron_lengths, trans_lengths = \
#                 Counter(), Counter(), Counter()
#             for chrom, strand, exons in trans.values():
#                 tran_len = 0
#                 for i, exon in enumerate(exons):
#                     exon_len = exon[1]-exon[0]+1
#                     exon_lengths[exon_len] += 1
#                     tran_len += exon_len
#                     if i == 0:
#                         continue
#                     intron_lengths[exon[0] - exons[i-1][1]] += 1
#                 trans_lengths[tran_len] += 1

#             print('genes: {}, genes with multiple isoforms: {}'.format(
#                   len(genes), sum(len(v) > 1 for v in genes.values())),
#                   file=stderr)
#             print('transcripts: {}, transcript avg. length: {:d}'.format(
#                     len(trans), sum(trans_lengths.elements())/len(trans)),
#                   file=stderr)
#             print('exons: {}, exon avg. length: {:d}'.format(
#                     sum(exon_lengths.values()),
#                     sum(exon_lengths.elements())/sum(exon_lengths.values())),
#                   file=stderr)
#             print('introns: {}, intron avg. length: {:d}'.format(
#                     sum(intron_lengths.values()),
#                     sum(intron_lengths.elements())/sum(intron_lengths.values())),
#                   file=stderr)
#             print('average number of exons per transcript: {:d}'.format(
#                     sum(exon_lengths.values())/len(trans)),
#                   file=stderr)


#     def run(self):
#         """Main."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             for gff in gffs:
#                 gtf_file = self.workdir + "/" +\
#                     gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
#                 out_file = self.workdir + "/" +\
#                     gff.split("/")[-1].rsplit(".", 1)[0] + ".splice"
#                 self.extract_splice_sites(gtf_file, out_file)
#         else:
#             gff = os.path.abspath(self.gff_file)
#             gtf_file = self.workdir + "/" +\
#                 gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
#             out_file = self.workdir + "/" +\
#                 self.gff_file.split("/")[-1].rsplit(".", 1)[0] + ".splice"
#             self.extract_splice_sites(gtf_file, out_file)


# # @requires(Hisat)
# class RefNames(luigi.Task):
#     """Extract the name of chromosomes where reads were mapped in BAM file."""

#     map_ref = Parameter()

#     def output(self):
#         """Split Files."""
#         return luigi.LocalTarget(self.map_ref)

#     def run(self):
#         """Extract names of chromosome from SAM file."""
#         bam_file = self.outsam.split(".sam")[0] + ".bam"
#         sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
#         grep_cmd = "samtools view -h %s |\
#                     grep -v -e '^.SQ' |\
#                     grep -v -e 'PN:' |\
#                     grep -v -e 'SO:' |\
#                     awk -F '\t' '{print $3}' |\
#                     sort | uniq |\
#                     grep -v '*' > %s" % (sorted_bam_file, self.map_ref)
#         subprocess.Popen(grep_cmd, shell=True)


# # @inherits(HisatMapW)
# class GetRefNames(luigi.WrapperTask):
#     """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

#     def requires(self):
#         """A pipeline that runs from mapping to count for euk and prok."""
#         splist = [self.workdir + "/" +
#                   f for f in os.listdir(self.workdir) if f.endswith('.splice')]
#         if len(splist) > 1:
#             splice_file = ','.join(splist)
#         else:
#             splice_file = splist[0]
#         for samp, fastq in self.fastq_dic.items():
#             trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
#             map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
#             yield RefNames(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
#                                    trim_dir + "/" + samp + ".2.trimmed.fastq"],
#                            num_cpus=self.num_cpus,
#                            indexfile=self.indexfile,
#                            spliceFile=splice_file,
#                            outsam=map_dir + "/" + samp + ".sam",
#                            map_ref=map_dir + "/" + samp + ".reflist",
#                            ref_file=self.ref_file,
#                            sample=samp,
#                            qc_outdir=trim_dir,
#                            map_dir=map_dir)




#@requires(Hisat)
class Split2ProkEuk(luigi.Task):
    """ Split BAM file to prok and euk"""

    ref_file = Parameter()
    workdir = Parameter()
    outsam = Parameter()

    def output(self):
        """Split Files."""
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        return [luigi.LocalTarget(prok_out), luigi.LocalTarget(euk_out)]

    def split_aln_file(self, fasta, bam_file, out_bamfile):
        """Split"""
        bed_file = os.path.join(self.workdir, "processes", "novel",
                                os.path.basename(fasta) + ".bedfile")
        # bed_file = workdir + "/" + os.path.basename(fasta) + ".bedfile"
        samtools_opt = ["view", "-b", "-L", bed_file, bam_file]
        samtools_cmd = (samtools[samtools_opt]) > out_bamfile
        samtools_cmd()

    def run(self):
        """ Split BAM file to proks and euks"""
        fastas = self.ref_file.split(",")
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        self.split_aln_file(fastas[0], bam_file, prok_out)
        self.split_aln_file(fastas[1], bam_file, euk_out)


#@inherits(HisatMapW)
class Split2ProkEukW(luigi.WrapperTask):
    """Group chromosomes to prokaryotic and eukaryotic."""
    fastq_dic = DictParameter()
    ref_file = Parameter()
    workdir = Parameter()

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        for samp, fastq in self.fastq_dic.items():
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            yield Split2ProkEuk(outsam=map_dir + "/" + samp + ".sam",
                                ref_file=self.ref_file,
                                workdir=self.workdir)
