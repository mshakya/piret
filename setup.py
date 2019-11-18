#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="piret",
    version='0.1.0',
    description="RNAseq pipeline",
    author="Migun Shakya",
    author_email="migun@lanl.gov",
    url="https://github.com/mshakya/piret",
    install_requires=open("requirements.txt").read().splitlines(),
    packages=find_packages(),
    scripts=['bin/piret', 'scripts/EdgeR', "scripts/plot_pathway",
             "scripts/RDESeq2", "scripts/gage_analysis",
             "scripts/Rballgown",
             "thirdparty/omics-pathway-viewer/scripts/opaver.r"],
    data_files=[("thirdparty/eggnog-mapper/",
                 ["thirdparty/eggnog-mapper/emapper.py",]),(
                 "thirdparty/eggnog-mapper/eggnogmapper",
                 ["thirdparty/eggnog-mapper/eggnogmapper/common.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/annota.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/orthology.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/search.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/server.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/seqio.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/utils.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/vars.py",
                 "thirdparty/eggnog-mapper/eggnogmapper/version.py"]),
                 ("thirdparty/omics-pathway-viewer/scripts/", [
                  "thirdparty/omics-pathway-viewer/scripts/opaver.pl"
                 ])
                 ],
    license="Apache License 2.0",
    platforms="Posix; MacOS X",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Bioinformatics",
    ],
)
