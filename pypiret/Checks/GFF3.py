#! /usr/bin/env python

"""Check design."""

import sys
import re
import collections
import pandas as pd

class CheckGFF():
    """Check experimental design are tab delimited."""

    def __init__(self, gff_file):
        """Initialize."""
        self.gff_col = pd.read_table(gff_file, sep='\t', comment="#")
        self.gff_col.columns = ["seqid", "source", "type", "start", "end",
                                "score", "strand", "phase", "attributes"]
    def size(self):
        """Check if gff file has 9 columns."""
        if len(self.gff_col.columns) != 9:
            sys.exit("""There must be 9 columns in a GFF3 file \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        else:
            return True

    def check_id(self):
        """Check if every row has ID in the 9th column."""
        gffw_child = self.gff_col.loc[self.gff_col['type'].isin(["mRNA", "gene", "CDS", "exon"])]
        nrows = gffw_child.shape[0]
        id_count = gffw_child.iloc[:, 8].str.contains('ID=').sum()
        if id_count == nrows:
            return True
        else:
            sys.exit("""PiReT requires mRNA, gene, CDS, and exon to have ID in atrributes\n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")

    def check_unique_id(self):
        """Check if every IDs are unique."""

        gff_mrna = self.gff_col[self.gff_col['type'].isin(["mRNA"])]
        gff_cds = self.gff_col[self.gff_col['type'].isin(["CDS"])]
        gff_gene = self.gff_col[self.gff_col['type'].isin(["gene"])]
        gff_exon = self.gff_col[self.gff_col['type'].isin(["exon"])]
        p_reg = re.compile('ID=.*?;')

        if gff_mrna.empty is False:
            gff_mrna['mRNA_ID'] = gff_mrna.apply(lambda row: p_reg.search(row['attributes']).group(0),
                                                 axis=1)
            mrna_ids = gff_mrna['mRNA_ID'].tolist()
        if gff_cds.empty is False:
            gff_cds['cds_ID'] = gff_cds.apply(lambda row: p_reg.search(row['attributes']).group(0),
                                              axis=1)
            cds_ids = gff_cds['cds_ID'].tolist()
        if gff_gene.empty is False:
            gff_gene['gene_ID'] = gff_gene.apply(lambda row: p_reg.search(row['attributes']).group(0), axis=1)
            gene_ids = gff_gene['gene_ID'].tolist()
        if gff_exon.empty is False:
            gff_exon['exon_ID'] = gff_exon.apply(lambda row: p_reg.search(row['attributes']).group(0), axis=1)
            exon_ids = gff_exon['exon_ID'].tolist()
        if gff_mrna.empty is False and gff_cds.empty is False:
            if len(list(set(mrna_ids).intersection(cds_ids))) > 0:
                print(set(mrna_ids).intersection(cds_ids))
                sys.exit("""Found overlapping IDs in mRNA and CDS \n
                            see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        if gff_mrna.empty is False and gff_gene.empty is False:
            print(set(mrna_ids).intersection(gene_ids))
            if len(list(set(mrna_ids).intersection(gene_ids))) > 0:
                sys.exit("""Found overlapping IDs in mRNA and gene \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        if gff_mrna.empty is False and gff_exon.empty is False:
            print(set(mrna_ids).intersection(exon_ids))
            if len(list(set(mrna_ids).intersection(exon_ids))) > 0:
                sys.exit("""Found overlapping IDs in mRNA and gene \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        if gff_exon.empty is False and gff_gene.empty is False:
            if len(list(set(gene_ids).intersection(exon_ids))) > 0:
                print(set(exon_ids).intersection(gene_ids))
                sys.exit("""Found overlapping IDs in gene and exon \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        if gff_gene.empty is False and gff_cds.empty is False:
            if len(list(set(gene_ids).intersection(cds_ids))) > 0:
                print(set(cds_ids).intersection(gene_ids))
                sys.exit("""Found overlapping IDs in gene and CDS \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")
        if gff_exon.empty is False and gff_cds.empty is False:            
            if len(list(set(exon_ids).intersection(cds_ids))) > 0:
                print(set(exon_ids).intersection(cds_ids))
                sys.exit("""Found overlapping IDs in exon and CDS \n
                        see: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md""")

        return True
