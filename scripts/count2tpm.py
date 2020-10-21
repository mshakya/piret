import re
import pandas as pd
import numpy as np


def read_counts2tpm(df, sample_name):
    """
    convert read counts to TPM (transcripts per million)
    :param df: a dataFrame contains the result coming from featureCounts
    :param sample_name: a list, all sample names, same as the result of featureCounts
    :return: TPM
    """
    result = df
    sample_reads = result.loc[:, sample_name].copy()
    gene_len = result.loc[:, ['Length']]
    rate = sample_reads.values / gene_len.values
    tpm = rate / np.sum(rate, axis=0).reshape(1, -1) * 1e6
    out_df = pd.DataFrame(data=tpm, columns=sample_name)
    out_df['Geneid'] = df['Geneid']
    return out_df


a = pd.DataFrame(data={
    'Gene': ("A", "B", "C", "D", "E"),
    'Length': (100, 50, 25, 5, 1),
    'S1': (80, 10,  6,  3,   1),
    'S2': (20, 20, 10, 50, 400)
})

read_counts2tpm(a, ['S1', 'S2'])
a.head()
a = pd.read_csv(
    "../covid_drugs/processes/featureCounts/prokarya/gene_count.tsv", sep="\t", skiprows=1)
a.columns = list(map(lambda x: re.sub(".*/", "", x), a.columns))

read_counts2tpm(a, ['2102_001_srt.bam'])
read_counts2tpm(a, ['2102_009_srt.bam'])
read_counts2tpm(a, ['2102_009_srt.bam', '2102_009_srt.bam'])
samp_list = [if '2102' in x for x in a.columns]
samp_list = [x if '2102' in x for x in a.columns]
samp_list = [x if '2102' in x for x in a.columns]
samp_list = [x for x in a.columns if '2102' in x]
samp_list
read_counts2tpm(a, samp_list)
%hist

read_counts2cpm(vir, samp_list)
mirna = pd.read_csv(
    "../covid_drugs/processes/featureCounts/eukarya/miRNA_count.tsv", sep="\t", skiprows=1)
mirna
mirna = pd.read_csv(
    "../covid_drugs/processes/featureCounts/eukarya/miRNA_name_count.tsv", sep="\t", skiprows=1)
mirna.columns = list(map(lambda x: re.sub(".*/", "", x), mirna.columns))
mirna
read_counts2tpm(mirna, samp_list)
mirna
mirna.columns = list(map(lambda x: re.sub(".*/", "", x), mirna.columns))
mirna
read_counts2tpm(mirna, samp_list)
mirna_tmp = read_counts2tpm(mirna, samp_list)
mirna_tmp.to_csv("mirna_tpm.csv", sep="\t", index=False)
pwd
