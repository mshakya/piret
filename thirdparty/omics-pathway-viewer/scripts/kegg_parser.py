#!/usr/bin/env python3
import numpy as np
import pandas as pd
import logging
import sqlite3
import re
import sys
import os

class KeggParser(object):
    """
        parsing kegg info from downloaded file
    """
    def __init__(self, db="kegg_info.db", log_level=logging.ERROR):
        self.logger = logging.getLogger('kegg_parser.kegg_parser')
        self.logger.setLevel(log_level)
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        
        self.db = db
        self.conn = None
        self.c = None
    
    def init_db(self, overwrite=False):
        if os.path.exists(self.db) and overwrite:
            os.remove(self.db)
        
        try:
            self.conn = sqlite3.connect(self.db)
            self.c = self.conn.cursor()

            self.c.execute("""
            CREATE TABLE IF NOT EXISTS kegg_pathway(
            ID         CHAR(5)     NOT NULL,
            NAME       CHAR(200)   NOT NULL,
            PRIMARY KEY (ID)
            );
            """)

            self.c.execute("""
            CREATE TABLE IF NOT EXISTS ec_info(
            ID        CHAR(20)    NOT NULL, /* ec id */
            NAME      CHAR(30)    NULL, /* ec name */
            CLASS     CHAR(30)    NULL, /* ec class */
            COMMENT   CHAR(200)   NULL,
            PATHWAY   CHAR(200)   NULL,
            PRIMARY KEY (ID)
            );
            """)

            self.c.execute("""
            CREATE TABLE IF NOT EXISTS ko_info(
            ID        CHAR(6)     NOT NULL, /* ko id */
            NAME      CHAR(100)   NULL, /* ko name */
            DEF       CHAR(200)   NULL, /* ko definition */
            PATHWAY   CHAR(200)   NULL,
            PRIMARY KEY (ID)
            );
            """)

            self.c.execute("""
            CREATE TABLE IF NOT EXISTS c_info(
            ID        CHAR(6)    NOT NULL, /* Compound id */
            NAME      CHAR(100)   NULL, /* Compound name */
            FORMULA   CHAR(100)   NULL,
            PATHWAY   CHAR(200)   NULL,
            PRIMARY KEY (ID)
            );
            """)
        except sqlite3.Error as e:
            self.logger.error("Database error: %s" % e)
        except Exception as e:
            self.logger.error("Exception in _query: %s" % e)

    def pwy_list_parser(self, info_file):
        with open(info_file) as f:
            for line in f:
                entry, description = line.rstrip().split("\t")
                # insert into database
                self.c.execute( 'INSERT INTO kegg_pathway VALUES ("%s", "%s")' % (entry.replace('path:map',''), description) )
                
    def info_parser(self, info_file):
        with open(info_file) as f:
            current_section = None
            entry = None
            name = None
            definition = None
            formula = None
            comment = None
            ec_class = None
            pwy_ids = []

            for line in f.readlines():
                # read through each lines
                section = line[:12].strip() # section names are within 12 columns
                if not section == "":
                    current_section = section

                if current_section == "NAME":
                    name = line[12:].strip()

                if current_section == "DEFINITION":
                    definition = line[12:].strip()

                if current_section == "CLASS":
                    ec_class = line[12:].strip()

                if current_section == "ENTRY":
                    entry = line[12:].replace("EC ","") # only apply to EC entry
                    entry = entry.split(" ")[0]
                    
                if current_section == "FORMULA":
                    formula = line[12:].strip()

                if current_section == "COMMENT":
                    comment = line[12:].strip()
                    
                if current_section == "PATHWAY":  # Get the pathways and add them to a list
                    pwy_fields = line[12:].split(" ")
                    pwy_ids.append( re.sub(r'^[a-z]+', '', pwy_fields[0]) )
            
            pathways = ", ".join(pwy_ids)

            self.logger.debug(f"file: {info_file}, entry={entry}, name={name}, ec_class={ec_class}, comment={comment}, definition={definition}, pathways={pathways}")

            if not entry:
                self.logger.warn(f"file: {info_file}, No entry found.")
                return

            if entry.startswith('C'):
                self.c.execute('INSERT INTO c_info VALUES (?, ?, ?, ?)', (entry, name, formula, pathways) )
            elif entry.startswith('E'):
                self.c.execute('INSERT INTO ec_info VALUES (?, ?, ?, ?, ?)', (entry, name, ec_class, comment, pathways) )
            elif entry.startswith('K'):
                self.c.execute('INSERT INTO ko_info VALUES (?, ?, ?, ?)', (entry, name, definition, pathways))
            
            f.close()

    def close(self):
        self.conn.commit()
        self.conn.close()
        self.logger.info(f"database committed and closed.")

if __name__ == "__main__":
    parser = KeggParser(log_level=logging.INFO)
    parser.init_db(overwrite=True)
    parser.pwy_list_parser("kegg_raw/kegg_pathway_list.txt")

    dir_path = "kegg_raw/kegg_info_ko/"
    for filename in os.listdir(dir_path):
        parser.info_parser(dir_path+filename)

    dir_path = "kegg_raw/kegg_info_ec/"
    for filename in os.listdir(dir_path):
        parser.info_parser(dir_path+filename)

    dir_path = "kegg_raw/kegg_info_cpd/"
    for filename in os.listdir(dir_path):
        parser.info_parser(dir_path+filename)

    parser.close()
    