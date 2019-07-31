#! /usr/bin/env python

"""parse config files."""
import os
import configparser


def check_aligner(aligner, hisat_index, workdir, star_db):
    """checks if the database for the aligner exist."""
    if aligner in ["hisat2", "hisat"]:
        if hisat_index is None:
            hisat_index_name = os.path.join(workdir, "processes", hisat_index)
        else:
            if os.path.exists(hisat_index) is False:
                hisat_index_name = os.path.join(workdir, "processes", "hisat_index")
            else:
                hisat_index_name = hisat_index
        return hisat_index_name
    elif aligner in ["star", "STAR"]:
        if star_db is None:
            new_stardb = os.path.join(workdir, "processes", "stardb")
        else:
            if os.path.exists(os.path.join(star_db, "chrName.txt")) is False:
                new_stardb = os.path.join(workdir, "processes", "stardb")
            else:
                new_stardb = star_db
        return new_stardb
