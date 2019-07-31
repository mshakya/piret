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
from luigi.interface import build
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


class RunEmapper(luigi.Task):
    """ Run emapper. """
    gff_file = luigi.Parameter()
    fasta_file = luigi.Parameter()
    workdir = luigi.Parameter()

    def output(self):
        """Expected output JSON."""
        jfile = os.path.join(self.workdir, "emapper", "emapper.txt")
        return luigi.LocalTarget(jfile)

    def run(self):
        """ Create fasta file."""
        self.gff2faa(self.gff_file, self.fasta_file)

    def gff2faa(self, gff_file, fasta):
        """reads in gff file and fasta to output proteome."""
        if os.path.exists(os.path.join(self.workdir, "databases")) is False:
            os.makedirs(os.path.join(self.workdir, "databases"))
        db_out = os.path.join(self.workdir, "databases", "piret.db")
        with open(os.path.join(self.workdir, "databases", "aas.faa"), "w") as f:
            if os.path.exists(db_out) is False:
                db = gffutils.create_db(gff_file, dbfn=db_out, force=True,
                                    keep_order=True,
                                    merge_strategy="create_unique")
            else:
                db = gffutils.FeatureDB(db_out, keep_order=True)
            for feat_obj in db.features_of_type("CDS"):
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







def gff2json(gff_file, fasta_file, out_json, translate=None):
    """A function that converts a gff file to JSON file."""

    # read in the gff file to a database
    db = gffutils.create_db(gff_file, dbfn="piret.db", force=True,
                            keep_order=True, merge_strategy="create_unique")

    feat_dic = {}  # an empty dictionary to append features
    with open(out_json, "w") as json_file:
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


# gff2faa("tests/data/GCF_000009065.1_ASM906v1_genomic.gff",
# "tests/data/GCF_000009065.1_ASM906v1_genomic.fna")
# gff2json("tests/data/GCF_000009065.1_ASM906v1_genomic.gff",
# "tests/data/GCF_000009065.1_ASM906v1_genomic.fna",
# "test.json", translate=translate)

# def count2json(jf, ):
#     """Reads in a JSON file and append new information like FPKMS, reads,
#     etc."""

#     with open(jf) as json_file:
#         data = json.load(json_file)
#         for p in data['people']:

