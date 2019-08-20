#! /usr/bin/env python

"""Check design."""
import gffutils
import simplejson
import json
from Bio.Seq import Seq
from Bio import SeqIO
import Bio
import os
import luigi
import pandas as pd
from luigi.interface import build
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from plumbum.cmd import time
from luigi.util import requires
from piret.Runs import FaQC


class GetAAs(luigi.Task):
    """Get amino acid sequences."""
    gff_file = luigi.Parameter()
    fasta_file = luigi.Parameter()
    workdir = luigi.Parameter()
    kingdom = luigi.Parameter()
    ave_map = luigi.FloatParameter()

    def requires(self):
        """Check if those two files are present."""
        for f in [self.gff_file, self.fasta_file]:
            return FaQC.RefFile(f)

    def gff2faa(self, gff_file, fasta):
        """reads in gff file and fasta to output proteome."""
        # get the list of CDS that past the threshold.
        imp_cds = self.get_imp_cds(self.ave_map)

        # make directories
        if os.path.exists(os.path.join(self.workdir, "processes", "databases")) is False:
            os.makedirs(os.path.join(self.workdir, "processes", "databases"))
        db_out = os.path.join(self.workdir, "processes", "databases", "piret.db")
        with open(os.path.join(self.workdir, "processes", "databases", "aas.faa"), "w") as f:
            if os.path.exists(db_out) is False:
                # create db if not already present
                db = gffutils.create_db(gff_file, dbfn=db_out, force=True,
                                        keep_order=True,
                                        merge_strategy="create_unique")
            else:
                # read db if its already present
                db = gffutils.FeatureDB(db_out, keep_order=True)
            for feat_obj in db.features_of_type("CDS"):
                if feat_obj.id in imp_cds.Geneid.tolist():
                    nt_seqs = feat_obj.sequence(self.fasta_file)
                    prot_seqs = self.translate(nt_seqs, "CDS")
                    try:
                        desc = feat_obj.attributes['product'][0]
                    except:
                        desc = "No annotation"
                    record = SeqRecord(Seq(prot_seqs, IUPAC.protein),
                                   id=feat_obj.id,
                                   description=desc)
                    SeqIO.write(record, f, "fasta")

    def translate(self, nucleotide, type):
        """Takes in a string of nucleotides and translate to AA."""
        if type == "CDS":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        elif type == "exon":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        else:
            aa = "not translated"
        return aa

    def get_imp_cds(self, ave_map):
        """Read in the read count table
            and filter out cds that do not make the threshold
            of average reads as shown in ave_map"""
        cds_table = os.path.join(self.workdir, "processes", "featureCounts",
                                 self.kingdom, "CDS_count.tsv")
        cds_df = pd.read_csv(cds_table, sep="\t", comment="#")

        cds_df["mean"] = cds_df.iloc[:, 6:].mean(axis=1)
        cds_df = cds_df[cds_df["mean"] > ave_map]
        return cds_df

    def run(self):
        """Create fasta file."""
        self.gff2faa(self.gff_file, self.fasta_file)

    def output(self):
        """Expected amino acid output."""
        aa_file = os.path.join(self.workdir, "processes", "databases",
                               "aas.faa")
        return luigi.LocalTarget(aa_file)


@requires(GetAAs)
class RunEmapper(luigi.ExternalTask):
    """ Run emapper.

        Get KEGG# and other"""
    query_coverage = luigi.FloatParameter()
    subject_coverage = luigi.FloatParameter()

    def run_emapper(self):
        """Using the amino acid fasta file, run emapper."""

        aa_file = os.path.join(self.workdir, "processes", "databases",
                               "aas.faa")
        egg_dir = os.path.join(self.workdir, "processes", "emapper", "emapper")
        if os.path.exists(egg_dir) is False:
            os.makedirs(egg_dir)

        emap = ["thirdparty/eggnog-mapper/emapper.py", "-i",
                aa_file, "-o", egg_dir, "--data_dir", "../eggnog-mapper/data/",
                "--dbtype", "seqdb", "-m", "diamond", "--target_orthologs", "one2one",
                "--query-cover", self.query_coverage,
                "--subject-cover", self.subject_coverage,
                "--temp_dir", egg_dir]
        time[emap]()

    def translate(self, nucleotide, type):
        """Takes in a string of nucleotides and translate to AA."""
        if type == "CDS":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        elif type == "exon":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        else:
            aa = "not translated"
        return aa

    def run(self):
        """ Create fasta file."""
        self.run_emapper()

    def output(self):
        """Expected output JSON."""
        jfile = os.path.join(self.workdir, "processes", "emapper",
                             "emapper.emapper.annotations")
        return luigi.LocalTarget(jfile)
