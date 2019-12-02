#! /usr/bin/env python

"""
First step of pipeline.

to check if files are correctly formatted
"""

import sys
import os
from piret import logger



# class Initialize(object):
#     """
#     First steps when starting the pipeline.

#     check
#     """

#     def __init__(self, workdir, kingdom, gff_prok, gff_euk,
#                  method):
#         """
#         Argument description.

#         __init__ arguments:
#             - ``workdir`` : working directory,
#             - ``kingdom: prokaryote or eukarya or both,
#             - ``gff_prok`` : prokaryoitic gff file
#             - ``gff_euk`` : eukaryotic gff file
#             - ``method`` : Deseq2 or edgeR
#         """
#         self.workdir = workdir
#         self.kingdom = kingdom
#         self.gff_prok = gff_prok
#         self.gff_euk = gff_euk
#         self.method = method

#     def create_dirs(self):
#         """
#         Create working directory.

#         checks if the permission
#         """
#         if not os.path.exists(self.workdir):
#             os.makedirs(self.workdir)
#             process = logger.create_logger(self.workdir, 'process', 'DEBUG')
#             os.makedirs('/'.join([self.workdir, 'sum_gene_count']))
#             os.makedirs('/'.join([self.workdir, 'tmp']))
#             os.makedirs('/'.join([self.workdir, 'differential_gene']))
#             if self.kingdom == 'both':
#                 os.makedirs('/'.join([self.workdir, 'differential_gene', 'eukarya']))
#                 os.makedirs('/'.join([self.workdir, 'differential_gene', 'prokaryote']))
#                 if self.method == 'both':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#                 elif self.method == 'DeSeq':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                 elif self.method == 'edgeR':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#             elif self.kingdom == 'prokaryote':
#                 os.makedirs('/'.join([self.workdir, 'differential_gene', 'prokaryote']))
#                 if self.method == 'both':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#                 elif self.method == 'DeSeq':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                 elif self.method == 'edgeR':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#             elif self.kingdom == 'eukarya':
#                 os.makedirs('/'.join([self.workdir, 'differential_gene', 'eukarya']))
#                 if self.method == 'both':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#                 elif self.method == 'DeSeq':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'DeSeq']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'DeSeq']))
#                 elif self.method == 'edgeR':
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'eukarya', 'edgeR']))
#                     os.makedirs('/'.join([self.workdir, 'differential_gene',
#                                           'prokaryote', 'edgeR']))
#             process.info("[Created directories]")
#         elif os.path.exists(self.workdir):
#             print("Directory already exists, Please remove the directory or name\
#     assign one that does not exist!")
#             sys.exit(1)
