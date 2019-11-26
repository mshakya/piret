#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.workflows import single_seq as ss
from plumbum.cmd import rm, cp

luigi.configuration.get_config().set('SummarizeQC', 'faqc_min_L', '50')
luigi.configuration.get_config().set('SummarizeQC', 'n_cutoff', '4')
luigi.configuration.get_config().set('Hisat', 'min_introlen', '20')
luigi.configuration.get_config().set('Hisat', 'max_introlen', '500000')
luigi.configuration.get_config().set('Hisat', 'rna_strandness', 'FR')
luigi.configuration.get_config().set('FeatureCounts', 'fid', 'ID')
luigi.configuration.get_config().set('FeatureCounts', 'stranded', '0')

class TestSS:

    def test_SingleSeq(self):
        """Test single seq."""
        # make it copy the luigi file first
        

        sseq = ss.SingleSeq(qc=True,
                            fastq_dic={"samp5":"tests/data/fastqs/BTT_test26_R1.fastq.gz:tests/data/fastqs/BTT_test26_R2.fastq.gz"},
                            aligner="hisat2",
                            ref_fasta="tests/data/test_prok.fna",
                            num_cpus=2,
                            local_scheduler=True,
                            hisat_index="tests/test_ss/h_index",
                            workdir="tests/test_ss",
                            kingdom="prokarya",
                            no_of_jobs=1,
                            p_value=0.05,
                            exp_desn_file="tests/exp_design_file.txt",
                            stardb_dir="star_db",
                            emap_dir="tests/",
                            gff_file=os.path.join("tests", "data", "test_prok.gff"),
                            pathway=True)
        sseq.run_faqc()
        sseq.create_db()
        sseq.map_reads({"samp5":"tests/test_ss/qc/samp5/samp5.1.trimmed.fastq:tests/test_ss/qc/samp5/samp5.1.trimmed.fastq"})
        sseq.map_summarize()
        sseq.extract_pp()
        sseq.NovelRegions()
        sseq.create_new_gff()
        gff=os.path.join("tests", "test_ss",
                                    "processes", "novel", "updated.gff")
        sseq.feature_count(new_gff=gff)
        sseq.merge_stringtie(new_gff=gff)
        assert os.stat(os.path.join("tests", "test_ss",
                                    "processes", "novel", "updated.gff")).st_size > 1
        rm_cmd1 = rm["-rf", "tests/test_ss"]
        rm_cmd1()

