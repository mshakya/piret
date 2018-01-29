#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import sys
import luigi
from luigi import LocalTarget
from luigi.contrib.external_program import ExternalProgramTask
from pypiret import Summ
from luigi.util import inherits, requires
from plumbum.cmd import Rscript


class RefFile(luigi.ExternalTask):
    """An ExternalTask like this."""

    path = luigi.Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


@requires(Summ.FeatureCounts)
class EdgeR(luigi.Task):
    """Find DGE using EdgeR."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using EdgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_file = edger_dir + "/" + file.split(".csv")[0] + "_RPKM" + ".csv"
                    return LocalTarget(out_file)

    def run(self):
        """Run EdgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for root, dirs, files in os.walk(fcount_dir):
            print(files)
            for file in files:
                if file.endswith("csv"):
                    name = file.split("_")[-2]
                    EdgeR_list =  [self.bindir + "/../scripts/EdgeR", "-r", os.path.join(root, file),
                              "-e", self.exp_design,
                              "-p", self.p_value,
                              "-n", name,
                              "-o", edger_dir]
                    EdgeR_cmd = Rscript[EdgeR_list]
                    EdgeR_cmd()

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': self.bindir + "/../scripts/" + ":" + os.environ["PATH"]}


# class EdgeR(ExternalProgramTask):
#     """Find DGE using EdgeR."""

#     exp_design = luigi.Parameter()
#     p_value = luigi.FloatParameter()

#     def requires(self):
#         """Require count file to be present."""
#         return RefFile(os.path.join(self.workdir, "featureCounts", "CDS.count"))

#     def output(self):
#         """Expected output of DGE using EdgeR."""
#         edger_dir = self.workdir + "/edgeR/" + self.kingdom
#         return LocalTarget(edger_dir + "/" + "CPM.csv")

#     def program_args(self):
#         """Run EdgeR."""
#         fcount_dir = self.workdir + "/featureCounts"
#         edger_dir = self.workdir + "/edgeR/" + self.kingdom
#         if not os.path.exists(edger_dir):
#             os.makedirs(edger_dir)
#         for root, dirs, files in os.walk(fcount_dir):
#             for file in files:
#                 if file.endswith("count"):
#                     return ["Rscript", self.bindir + "/../scripts/EdgeR",
#                             "-r", os.path.join(root, file),
#                             "-e", self.exp_design,
#                             "-p", self.p_value,
#                             "-o", edger_dir]

#     def program_environment(self):
#         """Environmental variables for this program."""
#         return {'PATH': self.bindir + "/../scripts/" + ":" + os.environ["PATH"]}
