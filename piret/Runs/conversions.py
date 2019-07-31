#! /usr/bin/env python

"""Check design."""
import os
import sys
import luigi
import shutil
from luigi import LocalTarget
from piret.Runs import Summ
from luigi.util import inherits, requires
import pandas as pd
import gffutils
DIR = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(DIR, "../../scripts"))
os.environ["PATH"] += ":" + script_dir
sys.path.insert(0, script_dir)
import logging
import simplejson


class conversions(luigi.Task):
    """Convert gene count, RPKM, fold change table to GeneID or locus tag
        and also to ones that have EC# or KO# when available."""
    gff_file = luigi.Parameter()
    gene_count_table = luigi.Parameter()
    gene_RPKM_table = luigi.Parameter()
    gene_CPM_table = luigi.Parameter()
    gene_fc_table = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        for root, dirs, files in os.walk(edger_dir):
            for file in files:
                if file.endswith("__sig.csv"):
                    out_folder = file.split(".csv")[0]
                    out_filepath = os.path.join(edger_dir, out_folder, "greater.csv")
                    return LocalTarget(out_filepath)
    


    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for file in os.listdir(fcount_dir):
            if file.endswith("tsv"):
                name  = file.split("_")[-2]
                edger_list = ["-r", os.path.join(fcount_dir, file),
                          "-e", self.exp_design,
                          "-p", self.p_value,
                          "-n", name,
                          "-o", edger_dir]
                #TODO: get the output that has locus tag
                edger_cmd = EdgeR[edger_list]
                logger = logging.getLogger('luigi-interface')
                logger.info(edger_cmd)
                edger_cmd()
                if file == "gene_count.tsv":
                #TODO:convert the first column to locus tag
                    if self.pathway is True:
                        path_list = ["-d", edger_dir,
                                     "-m", "edgeR", "-c", self.org_code] # get pathway information
                        path_cmd = plot_pathway[path_list]
                        logger.info(path_cmd)
                        path_cmd()
                    if self.GAGE is True:
                        gage_list = ["-d", edger_dir, "-m",
                                     "edgeR", "-c", self.org_code]
                        gage_cmd = gage_analysis[gage_list]
                        logger.info(gage_cmd)
                        gage_cmd()
        self.summ_summ()


class conver2json(luigi.Task):
    """ Summarizes and converts all the results to one big json file."""
    gff_file = luigi.Parameter()
    fasta_file = luigi.Parameter()
    pathway = luigi.Parameter()

    def output(self):
        """Expected output JSON."""
        jfile = os.path.join(self.workdir, "out.json")
        return LocalTarget(jfile)

    def run(self):
        """ Create JSON files."""
        gff2json(self.out_json)

    def gff2json(self, out_json):
        """A function that converts a gff file to JSON file."""

        # read in the gff file to a database
        db_path = os.paht.join(self.workdir, "processes", "database", "piret.db")
        db_path = gffutils.create_db(gff_file, dbfn="piret.db", force=True,
                                keep_order=True, merge_strategy="create_unique")

        feat_dic = {}  # an empty dictionary to append features
        with open(self.out_json, "w") as json_file:
            for feat_obj in db.all_features():
                feat_dic['seqid'] = feat_obj.seqid
                feat_dic['id'] = feat_obj.id
                feat_dic['source'] = feat_obj.source
                feat_type = feat_obj.featuretype
                feat_dic['featuretype'] = feat_type
                feat_dic['start'] = feat_obj.start
                feat_dic['end'] = feat_obj.end
                feat_dic['strand'] = feat_obj.strand
                feat_dic['frame'] = feat_obj.frame
                feat_dic['attributes'] = feat_obj.attributes.items()
                feat_dic['extra'] = feat_obj.extra
                if feat_dic['featuretype'] not in ["region"]:
                    nt_seqs = feat_obj.sequence(fasta_file)
                    feat_dic['nt_seq'] = nt_seqs
                    if feat_dic['featuretype'] in ["CDS"]:
                        feat_dic['aa_seqs'] = translate(nt_seqs, "CDS")
                feat_json = simplejson.dumps(feat_dic)
                simplejson.dump(feat_dic, json_file)
