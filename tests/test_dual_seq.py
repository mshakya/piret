#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.workflows import dual_seq as ds
from plumbum.cmd import rm, cp

luigi.configuration.get_config().set('SummarizeQC', 'faqc_min_L', '50')
luigi.configuration.get_config().set('SummarizeQC', 'n_cutoff', '4')
luigi.configuration.get_config().set('Hisat', 'min_introlen', '20')
luigi.configuration.get_config().set('Hisat', 'max_introlen', '500000')
luigi.configuration.get_config().set('Hisat', 'rna_strandness', 'FR')
luigi.configuration.get_config().set('FeatureCounts', 'fid', 'ID')
luigi.configuration.get_config().set('FeatureCounts', 'stranded', '0')
luigi.configuration.get_config().set('map_star', 'align_intron_min', '21')
luigi.configuration.get_config().set('map_star', 'align_intron_max', '0')

def test_DualSeq():
    """Test single seq."""
    # make it copy the luigi file first
    

    dseq = ds.DualSeq(qc=True,
                        fastq_dic={"samp5":"tests/data/fastqs/samp5.1.trimmed.fastq:tests/data/fastqs/samp5.2.trimmed.fastq",
                                    "samp6":"tests/data/fastqs/samp6.1.trimmed.fastq:tests/data/fastqs/samp6.2.trimmed.fastq"},
                        aligner="hisat2",
                        prok_fasta="tests/data/test_prok.fna",
                        prok_gff="tests/data/test_prok.gff",
                        euk_fasta="tests/data/eukarya_test.fa",
                        euk_gff="tests/data/eukarya_test.gff3",
                        num_cpus=2,
                        local_scheduler=True,
                        hisat_index="tests/test_ds/h_index",
                        workdir="tests/test_ds",
                        kingdom="both",
                        no_of_jobs=1,
                        p_value=0.05,
                        exp_desn_file="tests/exp_design_file.txt",
                        stardb_dir="star_db",
                        emap_dir="tests/",
                        gff_file=os.path.join("tests", "data", "test_prok.gff"),
                        pathway=True)
    dseq.run_faqc()
    dseq.create_db()
    dseq.map_reads()
    dseq.map_summarize()
    dseq.split_prokeuk()
    dseq.extract_pp()
    dseq.NovelRegions()
    dseq.create_new_gff()
    euk_gff=os.path.join("tests", "test_ds",
                                "processes", "novel", "euk_updated.gff")
    prok_gff=os.path.join("tests", "test_ds",
                                "processes", "novel", "prok_updated.gff")
    dseq.feature_counts(new_gff=euk_gff, kingdom="eukarya")
    dseq.feature_counts(new_gff=prok_gff, kingdom="prokarya")
    dseq.merge_stringtie(new_gff=prok_gff + "," + euk_gff)
    assert os.stat(os.path.join("tests", "test_ds",
                                "processes", "featureCounts", "eukarya",
                                 "CDS_count.tsv")).st_size > 1
    rm_cmd1 = rm["-rf", "tests/test_ds"]
    rm_cmd1()



def test_DualSeq_star():
    """Test single seq."""
    dseq = ds.DualSeq(qc=True,
                        fastq_dic={"samp5":"tests/data/BTT_test15_R1.1000.fastq:tests/data/BTT_test15_R2.1000.fastq"},
                        aligner="STAR",
                        prok_fasta="tests/data/test_prok.fna",
                        prok_gff="tests/data/test_prok.gff",
                        euk_fasta="tests/data/eukarya_test.fa",
                        euk_gff="tests/data/eukarya_test.gff3",
                        num_cpus=2,
                        local_scheduler=True,
                        hisat_index="tests/test_ds/h_index",
                        workdir="tests/test_ds_star",
                        kingdom="both",
                        no_of_jobs=1,
                        p_value=0.05,
                        exp_desn_file="tests/exp_design_file.txt",
                        stardb_dir="tests/test_ds_star/star_db",
                        emap_dir="tests/",
                        gff_file=os.path.join("tests", "data", "test_prok.gff"),
                        pathway=True)
    dseq.run_faqc()
    dseq.create_db()
    dseq.map_reads()
    dseq.map_summarize()
    dseq.split_prokeuk()
    dseq.extract_pp()
    dseq.NovelRegions()
    dseq.create_new_gff()
    euk_gff=os.path.join("tests", "test_ds_star",
                                "processes", "novel", "euk_updated.gff")
    prok_gff=os.path.join("tests", "test_ds_star",
                                "processes", "novel", "prok_updated.gff")
    dseq.feature_counts(new_gff=euk_gff, kingdom="eukarya")
    dseq.feature_counts(new_gff=prok_gff, kingdom="prokarya")
    assert os.stat(os.path.join("tests", "test_ds_star",
                                "processes", "featureCounts", "prokarya",
                                 "CDS_count.tsv")).st_size > 1
    rm_cmd1 = rm["-rf", "tests/test_ds_star"]
    rm_cmd1()