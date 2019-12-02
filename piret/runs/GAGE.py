# #! /usr/bin/env python

# """Check design."""
# import luigi
# from luigi import LocalTarget
# from piret import Summ
# from luigi.util import inherits, requires
# import pandas as pd
# from plumbum.cmd import Rscript


# @requires(DGE.edgeR)
# class GAGE(luigi.Task):
#     """Perform GAGE analysis."""

#     exp_design = luigi.Parameter()
#     p_value = luigi.FloatParameter()
#     bindir = luigi.Parameter()

#     def output(self):
#         """ Expected output of GAGE analysis."""

#     def run(self):
#         # if se
#         edger_dir =  os.path.join(self.workdir, "edgeR", self.kingdom)




#     def output(self):
#         """Expected output of DGE using edgeR."""
#         fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
#         edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
#         for root, dirs, files in os.walk(fcount_dir):
#             for file in files:
#                 if file.endswith("csv"):
#                     out_filename = file.split(".tsv")[0] + "_RPKM.csv"
#                     out_filepath = os.path.join(edger_dir, out_filename)
#                     return LocalTarget(out_filepath)

#     def run(self):
#         """Run edgeR."""
#         fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
#         edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
#         edger_location = os.path.join(self.bindir, "../scripts/edgeR.R")
#         if not os.path.exists(edger_dir):
#             os.makedirs(edger_dir)
#         for root, dirs, files in os.walk(fcount_dir):
#             for file in files:
#                 if file.endswith("tsv"):
#                     name = file.split("_")[-2]
#                     edger_list = [edger_location,
#                                   "-r", os.path.join(root, file),
#                                   "-e", self.exp_design,
#                                   "-p", self.p_value,
#                                   "-n", name,
#                                   "-o", edger_dir]
#                     edger_cmd = Rscript[edger_list]
#                     edger_cmd()
#         self.summ_summ()

#     def summ_summ(self):
#         """Summarize the summary table to be displayed in edge"""
#         edger_dir = self.workdir + "/edgeR/" + self.kingdom
#         all_files = os.listdir(edger_dir)
#         if all_files:
#             out_file = os.path.join(edger_dir, "summary_updown.csv")
#             summ_files = [pd.read_csv(os.path.join(edger_dir, file),
#                                   index_col=0) for file in all_files if "summary.csv" in file ]
#             summ_df = pd.concat(summ_files)
#             summ_df.to_csv(out_file)