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
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
from functools import reduce
from piret.miscs import RefFile

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
                edger_cmd = EdgeR[edger_list]
                logger = logging.getLogger('luigi-interface')
                logger.info(edger_cmd)
                edger_cmd()
                if file == "gene_count.tsv":
                    pass
        if self.pathway is True:
            path_list = ["-d", edger_dir,
                         "-m", "edgeR", "-c",
                         self.org_code]
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
    """Summarizes and converts all the results to one big JSON file."""
    gff_file = luigi.Parameter()
    fasta_file = luigi.Parameter()
    pathway = luigi.BoolParameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    method = luigi.ListParameter()
    NovelRegions = luigi.BoolParameter()

    def requires(self):
        flist = []
        if "edgeR" in self.method:
            cpm_file = os.path.join(self.workdir, "processes", "edgeR",
                                    self.kingdom, "gene",
                                    "gene_count_CPM.csv")
            flist.append(cpm_file)
        elif "DESeq2" in self.method:
            fpm_file = os.path.join(self.workdir, "processes", "DESeq2",
                                    self.kingdom, "gene",
                                    "gene_count_FPKM.csv")
            flist.append(fpm_file)
        return [RefFile(f) for f in flist]

    def output(self):
        """Expected output JSON."""
        if self.kingdom == "prokarya":
            jfile = os.path.join(self.workdir, "prokarya_out.json")
            return LocalTarget(jfile)
        elif self.kingdom == "eukarya":
            jfile = os.path.join(self.workdir, "eukarya_out.json")
            return LocalTarget(jfile)

    def run(self):
        """Create JSON files."""
        if self.kingdom == "prokarya":
            jfile = os.path.join(self.workdir, "prokarya_out.json")
        elif self.kingdom == "eukarya":
            jfile = os.path.join(self.workdir, "eukarya_out.json")
        self.gff2json(jfile)

    def gff2json(self, out_json):
        """A function that converts a gff file to JSON file."""
        db_dir = os.path.join(self.workdir, "processes", "databases", self.kingdom)
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)
        db_out = os.path.join(db_dir, "piret.db")
        if not os.path.exists(db_out):
            db = gffutils.create_db(self.gff_file, dbfn=db_out, force=True,
                                    keep_order=True,
                                    merge_strategy="create_unique")
        else:
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
        else:
            read_summ_NovelRegion = {}

        emaps = self.get_emapper()

        with open(out_json, "w") as json_file:
            json_file.write("[\n")
            first = True

            for feat_obj in db.all_features():
                feat_dic = {}
                feat_dic['seqid'] = feat_obj.seqid
                feat_dic['id'] = feat_obj.id
                feat_dic['source'] = feat_obj.source
                feat_type = feat_obj.featuretype
                feat_dic['featuretype'] = feat_type
                feat_dic['start'] = feat_obj.start
                feat_dic['end'] = feat_obj.end
                feat_dic['length'] = abs(feat_obj.end - feat_obj.start) + 1
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

                nt_obj = None
                if feat_type != "region":
                    try:
                        nt_seqs = feat_obj.sequence(self.fasta_file)
                        nt_obj = Seq(nt_seqs, generic_dna)
                        feat_dic['nt_seq'] = nt_seqs
                    except KeyError:
                        pass

                if feat_type == "CDS":
                    if nt_obj is not None:
                        feat_dic['aa_seqs'] = self.translate(nt_obj, "CDS")

                    self.assign_scores(feat_dic=feat_dic,
                                       edger_sdic=edger_summ_cds,
                                       deseq_sdic=deseq_summ_cds,
                                       feat_id=feat_obj.id)
                    feat_dic['read_count'] = read_summ_cds.get(feat_obj.id, None)

                    self.assign_dges(feat_type="CDS", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="edgeR", dge_dict=dge_edger_cds)
                    self.assign_dges(feat_type="CDS", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="DESeq2", dge_dict=dge_deseq_cds)

                    if emaps is not None:
                        feat_dic['emapper'] = emaps.get(feat_obj.id, None)
                    else:
                        feat_dic['emapper'] = None

                elif feat_type == "NovelRegion":
                    feat_dic['read_count'] = read_summ_NovelRegion.get(feat_obj.id, None)

                elif feat_type == "rRNA":
                    feat_dic['read_count'] = read_summ_rRNA.get(feat_obj.id, None)

                elif feat_type == "tRNA":
                    feat_dic['read_count'] = read_summ_tRNA.get(feat_obj.id, None)

                elif feat_type == "exon":
                    feat_dic['read_count'] = read_summ_exon.get(feat_obj.id, None)

                elif feat_type == "gene":
                    self.assign_scores(feat_dic=feat_dic,
                                       edger_sdic=edger_summ_genes,
                                       deseq_sdic=deseq_summ_genes,
                                       feat_id=feat_obj.id)
                    feat_dic['read_count'] = read_summ_gene.get(feat_obj.id, None)
                    feat_dic['ballgown_values'] = ballgown_gene_pm.get(feat_obj.id, None)
                    feat_dic['stringtie_values'] = stringtie_tpms.get(feat_obj.id, None)

                    self.assign_dges(feat_type="gene", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="edgeR", dge_dict=dge_edger_gene)
                    self.assign_dges(feat_type="gene", feat_dic=feat_dic,
                                     feat_id=feat_obj.id,
                                     method="DESeq2", dge_dict=dge_deseq_gene)

                if not first:
                    json_file.write(",\n")
                else:
                    first = False

                json.dump(feat_dic, json_file)

            json_file.write("\n]\n")

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
        """Get emapper result as a dictionary."""
        emapper_file = os.path.join(self.workdir, "processes", "emapper",
                                    self.kingdom,
                                    "emapper.emapper.annotations")
        if os.path.exists(emapper_file) is not True:
            return None

        emap = pd.read_csv(emapper_file, sep='\t', skiprows=[0, 1, 2],
                           skipinitialspace=True, skipfooter=3,
                           header=None, engine='python')

        # Set column names to the first row
        emap.columns = emap.iloc[0]
        # Drop the first row (since it's now the header) and set index
        emap = emap.drop(index=0).set_index('#query_name')

        result = emap.to_dict(orient="index")
        del emap
        return result

    def read_summary(self, feat_type):
        """Get read values as a dictionary."""
        read_file = os.path.join(self.workdir, "processes", "featureCounts",
                                 self.kingdom, feat_type + "_count_sorted.csv")
        if os.path.exists(read_file) is not True:
            return {}

        read_data = pd.read_csv(read_file, sep=",", index_col="Geneid")
        read_data.columns = [x.split(".mapping.")[-1].split(".")[0] for x in read_data.columns]
        cols_to_drop = [c for c in ["Unnamed: 0", "Chr", "Start", "End",
                                    "Strand", "Length", "total"]
                        if c in read_data.columns]
        read_data = read_data.drop(cols_to_drop, axis=1).astype(int)
        read_dict = read_data.to_dict(orient="index")
        del read_data
        return read_dict

    def dge_summary(self, feat_type, method):
        """Summarize SGE results from edgeR or DESeq2."""
        dge_dir = os.path.join(self.workdir, "processes", method,
                               self.kingdom, feat_type)
        dge_files = [f for f in glob.glob(dge_dir + "**/*et.csv", recursive=True)]
        dge_dicts = {}
        for file in dge_files:
            dge_df = pd.read_csv(file, sep=",", index_col=0)
            if method == "edgeR":
                cols_to_drop = [c for c in ["Geneid", "Chr", "Start", "End",
                                            "Strand", "Length"]
                                if c in dge_df.columns]
                dge_df = dge_df.drop(cols_to_drop, axis=1)
            dge_dicts[str(os.path.basename(file).replace(".csv", ""))] = dge_df.to_dict(orient="index")
            del dge_df
        return dge_dicts

    def assign_dges(self, feat_type, feat_dic, feat_id, method, dge_dict):
        """Assign dge values in JSON file."""
        dge_dir = os.path.join(self.workdir, "processes", method,
                               self.kingdom, feat_type)
        dge_files = [os.path.basename(f).replace(".csv", "")
                     for f in glob.glob(dge_dir + "**/*et.csv", recursive=True)]
        if len(dge_files) < 1:
            return
        for key in dge_dict:
            try:
                feat_dic[key + "__" + method] = dge_dict[key][feat_id]
            except KeyError:
                feat_dic[key + "__" + method] = None

    def pm_summary(self, feat_type, method):
        """Get CPM/FPM and RPKM/FPKM values."""
        if method == "edgeR":
            cpm_file = os.path.join(self.workdir, "processes", method,
                                    self.kingdom, feat_type,
                                    feat_type + "_count_CPM.csv")
            rpkm_file = os.path.join(self.workdir, "processes", method,
                                     self.kingdom, feat_type,
                                     feat_type + "_count_RPKM.csv")
        elif method == "DESeq2":
            cpm_file = os.path.join(self.workdir, "processes", method,
                                    self.kingdom, feat_type,
                                    feat_type + "_count_FPM.csv")
            rpkm_file = os.path.join(self.workdir, "processes", method,
                                     self.kingdom, feat_type,
                                     feat_type + "_count_FPKM.csv")
        else:
            return ({}, {})

        cpm_dict = {}
        rpkm_dict = {}

        if os.path.exists(cpm_file):
            cpm_df = pd.read_csv(cpm_file, sep=",", engine='python', index_col=0)
            cpm_dict = cpm_df.to_dict(orient="index")
            del cpm_df

        if os.path.exists(rpkm_file):
            rpkm_df = pd.read_csv(rpkm_file, sep=",", engine='python', index_col=0)
            rpkm_dict = rpkm_df.to_dict(orient="index")
            del rpkm_df

        return (cpm_dict, rpkm_dict)

    def pm_summary_ballgown(self):
        pm_file = os.path.join(self.workdir, "processes", "ballgown",
                               self.kingdom, "summpary_PMs.csv")
        if os.path.exists(pm_file) is not True:
            return {}
        pm_df = pd.read_csv(pm_file, sep=",", index_col=6)
        cols_to_drop = [c for c in ["t_id", "chr", "strand", "start", "end",
                                    "num_exons", "length", "gene_id", "gene_name"]
                        if c in pm_df.columns]
        pm_df = pm_df.drop(cols_to_drop, axis=1)
        pm_dict = pm_df.to_dict(orient="index")
        del pm_df
        return pm_dict

    def stringtie_tpm(self):
        """Get TPMs from stringtie."""
        stie_dir = os.path.join(self.workdir, "processes", "stringtie")
        stie_files = [f for f in glob.glob(stie_dir + "/**/*sTie.tab", recursive=True)]
        dflist = []
        for f in stie_files:
            df = pd.read_csv(f, sep="\t")
            cols_to_drop = [c for c in ["Gene Name", "Strand", "Start", "End"] if c in df.columns]
            df = df.drop(cols_to_drop, axis=1)
            samp_name = os.path.basename(f)
            samp = re.sub("_sTie.tab", "", samp_name)
            df.columns = ["GeneID", "Reference", samp + "_cov",
                          samp + "_FPKM", samp + "_TPM"]

            # CRITICAL FIX: Deduplicate BEFORE joining to prevent memory explosion
            df = df.drop_duplicates(subset=['GeneID', 'Reference'])
            # Set the index to allow efficient horizontal concatenation
            df = df.set_index(['GeneID', 'Reference'])

            dflist.append(df)

        if not dflist:
            return {}

        # Efficiently concatenate all dataframes horizontally instead of chaining merges
        finaldf = pd.concat(dflist, axis=1)
        del dflist

        finaldic = finaldf.to_dict(orient="index")
        del finaldf
        return finaldic

    def translate(self, nucleotide, type):
        """Takes in a string of nucleotides and translate to AA."""
        if type == "CDS":
            aa = nucleotide.translate()
        elif type == "exon":
            aa = nucleotide.translate()
        else:
            aa = "not translated"
        return str(aa)
