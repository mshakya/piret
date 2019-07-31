#! /usr/bin/env python

"""Check fasta."""
import Bio
import re
import pandas as pd
import sys


class CheckFasta():
    """Check different instances of fasta."""

    def __init__(self):
        """Initialize."""
        # self.design_file = design_file

    def confirm_fasta(self, fasta_file):
        """Check if the given file is a fasta."""
        with open(fasta_file, 'r') as f_in:
            # read the first line
            fasta_header = f_in.readline()
            # check if there is a >
            if fasta_header[:1] != ">":
                raise TypeError
            elif len(fasta_header.split(">")[1]) == 0:
                raise TypeError
            try:
                second_line = next(f_in)
                for x in list(set(second_line)):
                    if x in ['A', 'T', 'G', 'C', 'N']:
                        pass
                    else:
                        # TODO: check this for all sequences
                        raise TypeError
            except StopIteration:
                raise TypeError
            return True

    def duplicate(self, fasta_file):
        """Check if a fasta file has duplicate sequences (based on name)."""
        seqs = Bio.SeqIO.parse(fasta_file)
        seq_header = []
        for seq in seqs:
            seq_header.append(seq.description)
        if len(seq_header) != len(set(seq_header)):
            raise ValueError


def remove_special_chars(self, fasta_file, out_fasta):
    """Remove special characters from sequence name."""
    seqs = Bio.SeqIO.parse(fasta_file)
    for seq in seqs:
        seq_id = re.sub('[|]', '_', seq.id)
        seq.id = seq_id
        Bio.SeqIO.write(seq, out_fasta, "fasta")


def gff_fasta(fasta_file, gff_file):
    """Check if gff file and fasta file are compatible."""
    seqs = Bio.SeqIO.parse(fasta_file)
    df = pd.read_csv(gff_file, sep="\t", comment='#', header=None)

    seq_header = []
    for seq in seqs:
        fa_id = set(seq_header.append(seq.id))

    gff_id = set(df[0].tolist())
    if fa_id == gff_id:
        pass
    else:
        sys.exit('GFF file and Reference fasta do not have same ids!')


def match_gff_fasta(fasta_file, gff_file, out_fasta):
    """Rename the sequence baesd on the GFF file."""
    seqs = Bio.SeqIO.parse(fasta_file)
    df = pd.read_csv(gff_file, sep="\t", comment='#', header=None)

    acc_ids = set(df[0].tolist())

    for seq in seqs:
        for acc in acc_ids:
            if acc in seq.id:
                seq.description = acc + " " + seq.description
                Bio.SeqIO.write(seq, out_fasta, "fasta")


#TODO: implement these checks
def check_quotations(gff_file):
    """Check if there are quotations in gff file."""
    with open(gff_file) as g:
        content = g.readlines()
        for l in content:
            if '"' in l:
                print(l)
                raise TypeError("There is a quotation in above line, remove them first")
