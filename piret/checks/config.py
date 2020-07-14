#! /usr/bin/env python

"""parse config files."""
import os
import sys
import configparser


def parse_config(cfg_file):
    """
    parses config file.
    """

    # TODO: if a parameter is not present in config file, automatically assign a None to it
    # TODO: read boolean parameter as boolean, not as strings
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(cfg_file)
    para_dic = {}
    para_dic['stranded'] = config.get('DEFAULT', 'stranded')
    para_dic['aligner'] = config.get('DEFAULT', 'aligner')
    para_dic['kingdom'] = config.get('DEFAULT', 'kingdom')
    para_dic['threads'] = config.get('DEFAULT', 'threads')
    para_dic['fasta_euk'] = config.get('DEFAULT', 'fasta_euk')
    para_dic['gff_euk'] = config.get('DEFAULT', 'gff_euk')
    para_dic['fasta_prok'] = config.get('DEFAULT', 'fasta_prok')
    para_dic['gff_prok'] = config.get('DEFAULT', 'gff_prok')
    para_dic['emap_db'] = config.get('DEFAULT', 'emap_db')
    if config.get('DEFAULT', 'pathway') == "True":
        para_dic['pathway'] = True
    else:
        para_dic['pathway'] = False
    para_dic['method'] = config.get('DEFAULT', 'method')
    para_dic['q_value'] = config.get('DEFAULT', 'q_value')
    if config.get('DEFAULT', 'novel') == "True":
        para_dic['novel'] = True
    else:
        para_dic['novel'] = False
    para_dic['hisat_index'] = config.get('DEFAULT', 'hisat_index')
    para_dic['star_index'] = config.get('DEFAULT', 'star_index')
    if config.get('DEFAULT', 'qc') == "True":
        para_dic['qc'] = True
    return para_dic


def check_input(para_dic, parser):
    """ A function to check if complimentary options are selected. """
    # for prokarya
    if para_dic['kingdom'] in ["prokarya", "both"]:
        if para_dic["fasta_prok"] is None:
            parser.error("""prokarya requires fasta file to be specified using
                        fasta_prok in config file""")
        elif os.path.exists(para_dic["fasta_prok"]) is False:
            parser.error("""reference genome doesnt exist""")
        elif ',' in para_dic["fasta_prok"]:
            if ',' in para_dic["gff_prok"]:
                for fprok in para_dic["fasta_prok"].split(","):
                    if os.path.exists(fprok) is False:
                        parser.error("""reference genome doesnt exist""")
                for gprok in para_dic["gff_prok"].split(","):
                    if os.path.exists(gprok) is False:
                        parser.error("""reference genome doesnt exist""")
            else:
                if len(para_dic["fasta_prok"].split(",")) != len(para_dic["gff_prok"].split(",")):
                    parser.error("""Not all gff have """)
        if para_dic["gff_prok"] is None:
            parser.error("""prokarya requires gff file to be specified using
                        gff_prok in config file""")
        else:
            if os.path.exists(para_dic["gff_prok"]) is False:
                parser.error("""reference genome gff doesnt exist""")

    # for eukarya
    elif para_dic['kingdom'] in ["eukarya", "both"]:
        if para_dic["fasta_euk"] is None:
            parser.error("""eukarya requires fasta file to be specified using
                        fasta_euk in config file""")
        elif os.path.exists(para_dic["fasta_euk"]) is False:
            parser.error("""eukaryotic reference genome doesnt exist""")
        elif ',' in para_dic["fasta_euk"]:
            if ',' in para_dic["gff_euk"]:
                for fprok in para_dic["fasta_euk"].split(","):
                    if os.path.exists(fprok) is False:
                        parser.error("""reference genome doesnt exist""")
                for gprok in para_dic["gff_euk"].split(","):
                    if os.path.exists(gprok) is False:
                        parser.error("""reference genome doesnt exist""")
            else:
                if len(para_dic["fasta_euk"].split(",")) != len(para_dic["gff_euk"].split(",")):
                    parser.error("""Not all gff have """)
        if para_dic["gff_euk"] is None:
            parser.error("""eukarya requires gff file to be specified using
                        gff_euk in config file""")
        else:
            if os.path.exists(para_dic["gff_euk"]) is False:
                parser.error(
                    """eukaryotic reference genome gff doesnt exist""")

    # if "ballgown" in para_dic['method'] and para_dic["kingdom"] == "prokarya":
    #     sys.exit("""Ballgown does not work for prokaryotic genomes,
    #                 pick method that does not have ballgown!""")


def check_ref(ref_fasta=None, ref_gff=None):
    """Check if references are present."""
    if ',' in ref_fasta:
        if all([os.path.exists(f) for f in ref_fasta.split(",")]) is True:
            pass
        else:
            sys.exit("One of the reference fasta do not exist!")
    if ',' in ref_gff:
        if all([os.path.exists(f) for f in ref_gff.split(",")]) is True:
            return True
        else:
            sys.exit("One of the reference gff do not exist!")
    else:
        if os.path.exists(ref_gff) and os.path.exists(ref_fasta):
            if check_gff(ref_gff) is True:
                return True
        else:
            if os.path.exists(ref_fasta) is False:
                exit_message = ' '.join(
                    ("Reference FASTA", ref_fasta, "does not exist!"))
                sys.exit(exit_message)
            elif os.path.exists(ref_gff) is False:
                exit_message = ' '.join(
                    ("Reference GFF", ref_gff, "does not exist!"))
