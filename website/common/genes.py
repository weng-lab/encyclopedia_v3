#!/usr/bin/env python

import os, sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from db_utils import getcursor

class LookupGenes:
    def __init__(self, DBCONN):
        self.DBCONN = DBCONN
        self.tableNames = {"mm10" : "genes_mm10",
                           "hg19" : "genes_hg19"}

    def lookup(self, assembly, gene):
        with getcursor(self.DBCONN, "lookup") as curs:
            curs.execute("""
SELECT chrom, chromStart, chromEnd FROM {table}
WHERE lower(gene) = lower(%(gene)s)
""".format(table=self.tableNames[assembly]),
                             {"gene" : gene})
            if (curs.rowcount > 0):
                return curs.fetchall()
            return None

    def fuzzy_lookup(self, assembly, gene):
        with getcursor(self.DBCONN, "lookup") as curs:
            curs.execute("""
SELECT gene FROM {table}
WHERE gene ~ lower(%(gene)s)
""".format(table=self.tableNames[assembly]),
                             {"gene" : gene})
            if (curs.rowcount > 0):
                return [x[0] for x in curs.fetchall()]
            return None
