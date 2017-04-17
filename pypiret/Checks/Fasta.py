#! /usr/bin/env python

"""Check fasta."""
from Bio import SeqIO


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
                raise TypeError('Specified file is not a fasta')
            elif len(fasta_header.split(">")[1]) == 0:
                raise TypeError('Specified file is not a fasta')
            try:
                second_line = next(f_in)
                for x in list(set(second_line)):
                    if x in ['A', 'T', 'G', 'C', 'N']:
                        pass
                    else:
                        # TODO: check this for all sequences
                        raise TypeError('Illegel character in the sequence')
            except StopIteration:
                raise TypeError('Specified file is not a fasta')
            return True

    def duplicate(self, fasta_file):
        """Check if a fasta file has duplicate sequences (based on name)."""
        seqs = SeqIO.parse(fasta_file)
        seq_header = []
        for seq in seqs:
            seq_header.append(seq.description)
        if len(seq_header) != len(set(seq_header)):
            raise ValueError("duplicate sequences found")

