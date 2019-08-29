#! /usr/bin/env python

"""Check design."""
import os
import sys
import luigi
import shutil
from luigi import LocalTarget
from luigi.util import inherits, requires
import pandas as pd
import gffutils
import glob
DIR = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(DIR, "../../scripts"))
os.environ["PATH"] += ":" + script_dir
sys.path.insert(0, script_dir)
import logging
import json
import Bio
import re
from functools import reduce


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
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        out_filepath = os.path.join(edger_dir, "summary_updown.csv")
        return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for file in os.listdir(fcount_dir):
            if file.endswith("tsv"):
                name = file.split("_")[-2]
                edger_list = ["-r", os.path.join(fcount_dir, file),
                              "-e", self.exp_design,
                              "-p", self.p_value,
                              "-n", name,
                              "-o", edger_dir]
                # TODO: get the output that has locus tag
                edger_cmd = EdgeR[edger_list]
                logger = logging.getLogger('luigi-interface')
                logger.info(edger_cmd)
                edger_cmd()
                if file == "gene_count.tsv":
                    # TODO:convert the first column to locus tag
                    if self.pathway is True:
                        path_list = ["-d", edger_dir,
                                     "-m", "edgeR", "-c",
                                     self.org_code]  # get pathway information
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
    """ Summarizes and converts all the results to one big JSON file."""
    gff_file = luigi.Parameter()
    fasta_file = luigi.Parameter()
    pathway = luigi.BoolParameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    method = luigi.ListParameter()
    NovelRegions = luigi.BoolParameter()

    def output(self):
        """Expected output JSON."""
        if self.kingdom == "prokarya":
            jfile = os.path.join(self.workdir, "prokarya_out.json")
            return LocalTarget(jfile)
        elif self.kingdom == "eukarya":
            jfile = os.path.join(self.workdir, "eukarya_out.json")
            return LocalTarget(jfile)

    def run(self):
        """ Create JSON files."""
        if self.kingdom == "prokarya":    
            jfile = os.path.join(self.workdir, "prokarya_out.json")
        elif self.kingdom == "eukarya":
            jfile = os.path.join(self.workdir, "eukarya_out.json")

        self.gff2json(jfile)

    def gff2json(self, out_json):
        """A function that converts a gff file to JSON file."""
        # read in the gff file to a database
        if os.path.exists(os.path.join(self.workdir, "processes",
                                       "databases")) is False:
            os.makedirs(os.path.join(self.workdir, "processes",
                                     "databases"))
        db_out = os.path.join(self.workdir, "processes", "databases",
                              self.kingdom,
                              "piret.db")
        if os.path.exists(db_out) is False:
            # create db if not already present
            db = gffutils.create_db(self.gff_file, dbfn=db_out, force=True,
                                    keep_order=True,
                                    merge_strategy="create_unique")
        else:
            # read db if its already present
            db = gffutils.FeatureDB(db_out, keep_order=True)

        if "edgeR" in self.method:
            edger_summ_cds = self.pm_summary("CDS", "edgeR")
            edger_summ_genes = self.pm_summary("gene", "edgeR")
            dge_edger_cds = self.dge_summary("CDS", "edgeR")
            dge_edger_gene = self.dge_summary("gene", "edgeR")
        else:
            edger_summ_cds = ({}, {})
            edger_summ_genes = ({}, {})
            dge_edger_cds = {}
            dge_edger_gene = {}
        if "DESeq2" in self.method:
            deseq_summ_cds = self.pm_summary("CDS", "DESeq2")
            deseq_summ_genes = self.pm_summary("gene", "DESeq2")
            dge_deseq_cds = self.dge_summary("CDS", "DESeq2")
            dge_deseq_gene = self.dge_summary("gene", "DESeq2")
        else:
            deseq_summ_cds = ({}, {})
            deseq_summ_genes = ({}, {})
            dge_deseq_cds = {}
            dge_deseq_gene = {}
        if "ballgown" in self.method:
            ballgown_gene_pm = self.pm_summary_ballgown()
        else:
            ballgown_gene_pm = {}
        stringtie_tpms = self.stringtie_tpm()
        read_summ_cds = self.read_summary("CDS")
        read_summ_gene = self.read_summary("gene")
        read_summ_rRNA = self.read_summary("rRNA")
        read_summ_tRNA = self.read_summary("tRNA")
        read_summ_exon = self.read_summary("exon")
        if self.NovelRegions is True:
            read_summ_NovelRegion = self.read_summary("NovelRegion")

        emaps = self.get_emapper()
        with open(out_json, "w") as json_file:
            json_list = []
            for feat_obj in db.all_features():
                feat_dic = {}  # an empty dictionary to append features
                feat_dic['seqid'] = feat_obj.seqid
                feat_dic['id'] = feat_obj.id
                feat_dic['source'] = feat_obj.source
                feat_type = feat_obj.featuretype
                feat_dic['featuretype'] = feat_type
                feat_dic['start'] = feat_obj.start
                feat_dic['end'] = feat_obj.end
                feat_dic["length"] = abs(feat_obj.end - feat_obj.start) + 1
                feat_dic['strand'] = feat_obj.strand
                feat_dic['frame'] = feat_obj.frame
                try:
                    feat_dic['locus_tag'] = feat_obj.attributes['locus_tag'][0]
                except KeyError:
                    pass
                try:
                    feat_dic['Note'] = feat_obj.attributes['Note']
                except KeyError:
                    pass
                feat_dic['extra'] = feat_obj.extra
                if feat_type != "region":
                    nt_seqs = feat_obj.sequence(self.fasta_file)
                    feat_dic['nt_seq'] = nt_seqs
# ============================================================================#
                if feat_type == "CDS":
                    # translate the CDS
                    feat_dic['aa_seqs'] = self.translate(nt_seqs, "CDS")

                    # assign FPKMs and FPMs
                    self.assign_scores(feat_dic=feat_dic,
                                       edger_sdic=edger_summ_cds,
                                       deseq_sdic=deseq_summ_cds,
                                       feat_id=feat_obj.id)
                    # assign read numbers
                    try:
                        feat_dic["read_count"] = read_summ_cds[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
                    # assign dge information
                    self.assign_dges(feat_type="CDS", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="edgeR", dge_dict=dge_edger_cds)
                    self.assign_dges(feat_type="CDS", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="DESeq2", dge_dict=dge_deseq_cds)
                    # assign EC#s, KOs, etc.
                    try:
                        feat_dic["emapper"] = emaps[feat_obj.id]
                    except KeyError:
                        feat_dic["emapper"] = None
# ============================================================================#
                elif feat_type == "NovelRegion":
                    try:
                        feat_dic["read_count"] = read_summ_NovelRegion[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
# ============================================================================#
                elif feat_type == 'rRNA':
                    try:
                        feat_dic["read_count"] = read_summ_rRNA[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
# ============================================================================#
                elif feat_type == 'tRNA':
                    try:
                        feat_dic["read_count"] = read_summ_tRNA[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
# ============================================================================#
                elif feat_type == 'exon':
                    try:
                        feat_dic["read_count"] = read_summ_exon[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
# ============================================================================#
                elif feat_type == "gene":
                    # assign scores
                    self.assign_scores(feat_dic=feat_dic,
                                       edger_sdic=edger_summ_genes,
                                       deseq_sdic=deseq_summ_genes,
                                       feat_id=feat_obj.id)
                    # assign read numbers
                    try:
                        feat_dic["read_count"] = read_summ_gene[feat_obj.id]
                    except KeyError:
                        feat_dic["read_count"] = None
                    # assign ballgown info
                    try:
                        feat_dic["ballgown_values"] = ballgown_gene_pm[feat_obj.id]
                    except KeyError:
                        feat_dic["ballgown_values"] = None
                    
                    # assign stringtie
                    try:
                        feat_dic["stringtie_values"] = stringtie_tpms[feat_obj.id]
                    except KeyError:
                        feat_dic["stringtie_values"] = None
                    # assign dge information
                    self.assign_dges(feat_type="gene", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="edgeR", dge_dict=dge_edger_gene)
                    self.assign_dges(feat_type="gene", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="DESeq2", dge_dict=dge_deseq_gene)
                else:
                    # print(feat_type)
                    pass
                # json.dump(feat_dic, json_file, indent=4)
                json_list.append(feat_dic)
            json.dump(json_list, json_file, indent=4)

    def assign_scores(self, feat_dic, edger_sdic, deseq_sdic, feat_id):
        """Assign scores from edger and deseq to summary dic."""
        try:
            feat_dic["edger_cpm"] = edger_sdic[0][feat_id]
        except KeyError:
            feat_dic["edger_cpm"] = None
        try:
            feat_dic["deseq_fpm"] = deseq_sdic[0][feat_id]
        except KeyError:
            feat_dic["deseq_fpm"] = None
        try:
            feat_dic["edger_rpkm"] = edger_sdic[1][feat_id]
        except KeyError:
            feat_dic["edger_rpkm"] = None
        try:
            feat_dic["deseq_fpkm"] = deseq_sdic[1][feat_id]
        except KeyError:
            feat_dic["deseq_fpkm"] = None

    def get_emapper(self):
        """get emapper result as a dataframe."""
        emapper_files = os.path.join(self.workdir, "processes", "emapper",
                                     self.kingdom,
                                     "emapper.emapper.annotations")
        if os.path.exists(emapper_files) is True:
            emap = pd.read_csv(emapper_files, sep='\t', skiprows=[0,1,2],
                               skipinitialspace=True, skipfooter=3,
                               header=None, engine='python')
            emap1 = emap.reset_index()
            emap1.columns = emap1.iloc[0]
            emap2 = emap1.drop(0).set_index('#query_name').to_dict(orient="index")
            return emap2
        else:
            return None

    def read_summary(self, feat_type):
        """Get read values as a dictionary."""
        read_file = os.path.join(self.workdir, "processes", "featureCounts",
                                 self.kingdom, feat_type + "_count_sorted.csv")

        if os.path.exists(read_file) is True:
            read_data = pd.read_csv(read_file, sep=",",
                                    index_col="Geneid")
            read_dict = read_data.drop(["Unnamed: 0", "Chr", "Start", "End",
                                        "Strand", "Length"], axis=1).to_dict(orient="index")
        else:
            read_dict = {}
        return read_dict

    def dge_summary(self, feat_type, method):
        """summarize SGE results from edgeR of DESeq2"""
        dge_dir = os.path.join(self.workdir, "processes", method,
                                self.kingdom, feat_type)
        dge_files = [f for f in glob.glob(dge_dir + "**/*et.csv", recursive=True)]
        dge_dicts = {}
        for file in dge_files:
            dge_df = pd.read_csv(file, sep=",", index_col=0)
            if method == "edgeR":
                dge_dict = dge_df.drop(["Geneid", "Chr", "Start", "End",
                                        "Strand", "Length"],
                                       axis=1).to_dict(orient="index")
            elif method == "DESeq2":
                dge_dict = dge_df.to_dict(orient="index")
            dge_dicts[str(os.path.basename(file).replace(".csv", ""))] = dge_dict
        return dge_dicts

    def assign_dges(self, feat_type, feat_dic, feat_id, method, dge_dict):
        """Assign dge values in JSON file."""
        dge_dir = os.path.join(self.workdir, "processes", method,
                               self.kingdom, feat_type)
        dge_files = [os.path.basename(f).replace(".csv", "")
                     for f in glob.glob(dge_dir + "**/*et.csv",
                     recursive=True)]
        if len(dge_files) < 1:
            pass
        else:
            for key, value in dge_dict.items():
                try:
                    feat_dic[key + "__" + method] = dge_dict[key][feat_id]
                except KeyError:
                    feat_dic[key + "__" + method] = None

    def pm_summary(self, feat_type, method):
        """Get FPKM values."""
        if method == "edgeR":
            cpm_file = os.path.join(self.workdir, "processes", method,
                                    self.kingdom,  feat_type,
                                    feat_type + "_count_CPM.csv")
            rpkm_file = os.path.join(self.workdir, "processes", method,
                                     self.kingdom,  feat_type,
                                     feat_type + "_count_RPKM.csv")
        elif method == "DESeq2":
            cpm_file = os.path.join(self.workdir, "processes", method,
                                    self.kingdom,  feat_type,
                                    feat_type + "_count_FPM.csv")
            rpkm_file = os.path.join(self.workdir, "processes", method,
                                     self.kingdom,  feat_type,
                                     feat_type + "_count_FPKM.csv")
        if all([os.path.exists(cpm_file), os.path.exists(rpkm_file)]) is False:
            return ({}, {})
        elif all([os.path.exists(cpm_file), os.path.exists(rpkm_file)]) is True:
            cpm_dict = pd.read_csv(cpm_file, sep=",",
                                   index_col=0).to_dict(orient="index")
            rpkm_dict = pd.read_csv(rpkm_file, sep=",",
                                    index_col=0).to_dict(orient="index")          
            return(cpm_dict, rpkm_dict)
        elif os.path.exists(cpm_file) is True and os.path.exists(rpkm_file) is False:
            cpm_dict = pd.read_csv(cpm_file, sep=",",
                                   index_col=0).to_dict(orient="index")
            return(cpm_dict, {})
        elif os.path.exists(rpkm_file) is True and os.path.exists(rpkm_file) is False:
            rpkm_dict = pd.read_csv(rpkm_file, sep=",",
                                    index_col=0).to_dict(orient="index")
            return({}, rpkm_dict)

    def pm_summary_ballgown(self):
        pm_file = os.path.join(self.workdir, "processes", "ballgown",
                               self.kingdom, "summpary_PMs.csv")
        if os.path.exists(pm_file) is True:
            pm_dict = pd.read_csv(pm_file, sep=",",
                                  index_col=6).drop(["t_id", "chr", "strand",
                                                      "start", "end",
                                                      "num_exons", "length",
                                                      "gene_id", "gene_name"],
                                                      axis=1).to_dict(orient="index")
            return pm_dict
        else:
            return {}

    def stringtie_tpm(self):
        """get TPMs from stringtie."""
        stie_dir = os.path.join(self.workdir, "processes", "stringtie")
        stie_files = [f for f in glob.glob(stie_dir + "/**/*sTie.tab",
                      recursive=True)]
        dflist = []
        for f in stie_files:
            df = pd.read_csv(f, sep="\t").drop(["Gene Name", "Reference", "Strand",
                                               "Start", "End"], axis=1)
            
            samp_name = os.path.basename(f)
            samp = re.sub("_.*", "", samp_name)
            df.columns = ["GeneID", samp + "_cov",
                          samp + "_FPKM", samp + "_TPM"]
            dflist.append(df)
            
        finaldf = reduce(lambda df1, df2: pd.merge(df1, df2, on='GeneID'), dflist)    
        finaldic = finaldf.set_index('GeneID').to_dict(orient="index")
        return finaldic

    def translate(self, nucleotide, type):
        """Takes in a string of nucleotides and translate to AA."""
        if type == "CDS":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        elif type == "exon":
            aa = Bio.Seq.translate(nucleotide, cds=False)
        else:
            aa = "not translated"
        return aa
