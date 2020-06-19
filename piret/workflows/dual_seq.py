import os
import sys
import argparse
import luigi
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
bin_path = os.path.join(lib_path, 'bin')
sys.path.append(lib_path)
os.environ["PATH"] += os.pathsep + bin_path
from piret.checks.design import CheckDesign
from piret.maps import hisat2
from piret.qc import FaQC
from piret.runs import Map
from piret.novel import srna
from piret.maps import star
from piret.counts import featurecounts
from piret.dge import edgeR
from piret.dge import DESeq2
from piret.dge import ballgown
from piret.counts import stringtie
from piret.summary import summarize
from piret.functions import function
from piret.pathways import opaver
from luigi.interface import build



class DualSeq:
    """Class that pieces luigi task for dual seq."""

    def __init__(self, **kwargs):
        self.qc = kwargs["qc"]
        self.fastq_dic = kwargs["fastq_dic"]
        self.workdir = kwargs["workdir"]
        self.num_cpus = kwargs["num_cpus"]
        self.aligner = kwargs['aligner']
        self.prok_fasta = kwargs['prok_fasta']
        self.euk_fasta = kwargs['euk_fasta']
        self.prok_gff = kwargs['prok_gff']
        self.euk_gff = kwargs['euk_gff']
        self.exp_desn_file = kwargs['exp_desn_file']
        self.pathway = kwargs["pathway"]
        self.hisat_index = kwargs["hisat_index"]        
        self.kingdom = kwargs["kingdom"]
        self.local_scheduler = kwargs["local_scheduler"]
        self.p_value = kwargs["p_value"]
        self.stardb_dir = kwargs["stardb_dir"]
        self.emap_dir = kwargs["emap_dir"]

    def run_faqc(self):
        """A function that calls the FaQC function.

        it returns QCed files in respective directory
        """
        if self.qc is True:
            build([FaQC.SummarizeQC(fastq_dic=self.fastq_dic,
                                    num_cpus=self.num_cpus,
                                    workdir=self.workdir)],
                  local_scheduler=self.local_scheduler, workers=1)
            qc_dic = {}
            for samp, path in self.fastq_dic.items():
                trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
                qc_dic[samp] = trim_dir + "/" + samp + ".1.trimmed.fastq" + ":" + \
                               trim_dir + "/" + samp + ".2.trimmed.fastq"          
            return qc_dic

    def create_db(self):
        """Function to create hisat index."""
        if self.aligner == "hisat2":
            ref_fastas = self.prok_fasta + "," + self.euk_fasta
            build([hisat2.HisatIndex(fasta=ref_fastas,
                                     hi_index=self.hisat_index,
                                     num_cpus=self.num_cpus),
                   Map.SAMindex(fasta=ref_fastas, workdir=self.workdir)],
                #    Map.CreateSplice(gff_file=self.ref_gffs.split(",")[1],
                            # workdir=self.workdir),
                #    Map.GetChromName(prok_ref=self.ref_fastas.split(",")[0],
                            # euk_ref=self.ref_fastas.split(",")[1],
                            # workdir=self.workdir)],
                          local_scheduler=self.local_scheduler,
                          workers=1)
        elif self.aligner in ["STAR", "star"]:
            ref_fastas = self.prok_fasta + "," + self.euk_fasta
            build([star.STARindex(fasta=ref_fastas,
                                num_cpus=self.num_cpus,
                                gff_file=self.euk_gff,
                                stardb_dir = self.stardb_dir,
                                kingdom=self.kingdom),
                   Map.SAMindex(fasta=ref_fastas, workdir=self.workdir)
                #    Map.CreateSplice(gff_file=self.euk_gff,
                #             workdir=self.workdir),
                #    Map.GetChromName(prok_ref=self.prok_fasta,
                #                     euk_ref=self.euk_fasta,
                #             workdir=self.workdir)
                            ],
                            local_scheduler=self.local_scheduler)

    def map_reads(self):
        """Function to map reads."""
        if self.aligner == "hisat2":
            build([hisat2.HisatMapW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index, workdir=self.workdir,
                                 kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler)
        elif self.aligner in ["STAR", "star"]:
            build([star.map_starW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 stardb_dir=self.stardb_dir, workdir=self.workdir)],
                  local_scheduler=self.local_scheduler)

    def map_summarize(self):
        """Summarize mapped reads into a table."""
        if self.aligner == "hisat2":
            build([hisat2.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    indexfile=self.hisat_index,
                                    num_cpus=self.num_cpus,
                                    kingdom=self.kingdom)],
            local_scheduler=self.local_scheduler, workers=1)
        elif self.aligner == "STAR":
            build([star.SummarizeStarMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    stardb_dir=self.stardb_dir,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)

    def split_prokeuk(self):
        build([Map.Split2ProkEukW(fastq_dic=self.fastq_dic,
                                workdir=self.workdir,
                                  ref_file=self.prok_fasta + "," + self.euk_fasta)],
                          local_scheduler=self.local_scheduler,
                          workers=1)

    def summarize_map(self):
        """Summarize mapped reads into a table."""
        if self.aligner == "hisat2":
            build([Map.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    indexfile=self.hisat_index,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)
        elif self.aligner == "STAR":
            build([Map.SummarizeStarMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    stardb_dir=self.stardb_dir,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)

    def extract_pp(self):
        """Extract properly paired reads."""
        build([srna.ExtractPPW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index, workdir=self.workdir,
                                 kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler)

    def NovelRegions(self):
        """Find novel regions."""
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                  workdir=self.workdir,
                                 kingdom=self.kingdom,
                                 gff_file=self.prok_gff + "," + self.euk_gff)],
              local_scheduler=self.local_scheduler)

    def create_new_gff(self):
        build([srna.CompileGFF(fastq_dic=self.fastq_dic,
                               kingdom=self.kingdom,
                               workdir=self.workdir,
                               gff_file=self.prok_gff + "," + self.euk_gff)],
          local_scheduler=self.local_scheduler, workers=1)

    def merge_stringtie(self, new_gff):
        build([stringtie.MergeStringTies(fastq_dic=self.fastq_dic,
                                         num_cpus=self.num_cpus,
                                workdir=self.workdir,
                                gff_file=new_gff,
                                kingdom=self.kingdom)],
           local_scheduler=self.local_scheduler, workers=1)
    
    def restringtie(self):
        build([stringtie.ReStringTieScoresW(fastq_dic=self.fastq_dic,
                                       num_cpus=self.num_cpus,
                                       workdir=self.workdir,
                                       kingdom=self.kingdom)],
           local_scheduler=self.local_scheduler, workers=1)

    def run_ballgown(self):
        build([ballgown.ballgown(
                            workdir=self.workdir,
                            kingdom="eukarya",
                            exp_design=self.exp_desn_file,
                            p_value=self.p_value)],
               local_scheduler=self.local_scheduler, workers=1)

    def feature_counts(self, new_gff, kingdom):
        build([featurecounts.FeatureCounts(fastq_dic=self.fastq_dic,
                                           num_cpus=self.num_cpus,
                                           gff_file=new_gff,
                                           indexfile=self.hisat_index,
                                           kingdom=kingdom,
                                           workdir=self.workdir,
                                           ref_file=self.prok_fasta + "," + self.euk_fasta)],
          local_scheduler=self.local_scheduler, workers=1)

    def run_edger(self, kingdom):
        build([edgeR.edgeR(workdir=self.workdir,
                           kingdom=kingdom,
                           exp_design=self.exp_desn_file,
                           p_value=self.p_value)], local_scheduler=self.local_scheduler,
                           workers=1)

    def run_deseq2(self, kingdom):
        build([DESeq2.DESeq2(workdir=self.workdir,
                             kingdom=kingdom,
                             exp_design=self.exp_desn_file,
                             p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=1)

    def run_emapper(self, new_gff, kingdom, fasta):
        build([function.RunEmapper(workdir=self.workdir,
                                   gff_file=new_gff,
                                   fasta_file=fasta,
                                   kingdom=kingdom,
                                   emapper_dir=self.emap_dir)],
              local_scheduler=self.local_scheduler, workers=1)

    def run_opaver(self, method, kingdom):
        build([opaver.RunOpaver(workdir=self.workdir,
                         kingdom=kingdom,
                         method=method)],
              local_scheduler=self.local_scheduler, workers=1)

    def summ_json(self, new_gff, method, NovelRegions, kingdom, fasta):
        build([summarize.conver2json(gff_file=new_gff,
                                     fasta_file=fasta,
                                     pathway=self.pathway,
                                     workdir=self.workdir,
                                     kingdom=kingdom,
                                     method=method,
                                     NovelRegions=NovelRegions)],
              local_scheduler=self.local_scheduler, workers=1)