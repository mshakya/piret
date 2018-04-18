#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi import LocalTarget
from pypiret import Summ
from luigi.util import inherits, requires
from plumbum.cmd import Rscript



@requires(Summ.FeatureCounts)
class edgeR(luigi.Task):
    """Find DGE using edgeR."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_file = edger_dir + "/" + file.split(".csv")[0] + "_RPKM.csv"
                    return LocalTarget(out_file)

    def run(self):
        """Run edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("tsv"):
                    name = file.split("_")[-2]
                    edger_list = [self.bindir + "/../scripts/edgeR.R",
                                  "-r", os.path.join(root, file),
                                  "-e", self.exp_design,
                                  "-p", self.p_value,
                                  "-n", name,
                                  "-o", edger_dir]
                    edger_cmd = Rscript[edger_list]
                    edger_cmd()

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': self.bindir + "/../scripts/" + ":" + os.environ["PATH"]}


# class edgeR(ExternalProgramTask):
#     """Find DGE using edgeR."""

#     exp_design = luigi.Parameter()
#     p_value = luigi.FloatParameter()

#     def requires(self):
#         """Require count file to be present."""
#         return RefFile(os.path.join(self.workdir, "featureCounts", "CDS.count"))

#     def output(self):
#         """Expected output of DGE using edgeR."""
#         edger_dir = self.workdir + "/edgeR/" + self.kingdom
#         return LocalTarget(edger_dir + "/" + "CPM.csv")

#     def program_args(self):
#         """Run edgeR."""
#         fcount_dir = self.workdir + "/featureCounts"
#         edger_dir = self.workdir + "/edgeR/" + self.kingdom
#         if not os.path.exists(edger_dir):
#             os.makedirs(edger_dir)
#         for root, dirs, files in os.walk(fcount_dir):
#             for file in files:
#                 if file.endswith("count"):
#                     return ["Rscript", self.bindir + "/../scripts/edgeR",
#                             "-r", os.path.join(root, file),
#                             "-e", self.exp_design,
#                             "-p", self.p_value,
#                             "-o", edger_dir]

#     def program_environment(self):
#         """Environmental variables for this program."""
#         return {'PATH': self.bindir + "/../scripts/" + ":" + os.environ["PATH"]}
