#!/usr/bin/env python

import os, sys, json, psycopg2, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils
from dbs import DBS
from files_and_paths import Dirs
from coord import Coord
from db_utils import getcursor

class Genes:
    def __init__(self, DBCONN):
        self.DBCONN = DBCONN
        self.tableNames = {"mm9" : "genes_mm9",
                           "mm10" : "genes_mm10",
                           "hg19" : "genes_hg19"}

    def lookup(self, assembly, gene):
        with getcursor(self.DBCONN, "lookup") as curs:
            curs.execute("""
SELECT chrom, chromStart, chromEnd FROM {table}
WHERE gene = %(gene)s
""".format(table=self.tableNames[assembly]),
                             {"gene" : gene})
            if (curs.rowcount > 0):
                return curs.fetchall()
            return None

def setupAndCopy(cur, fnp, table_name):
    print "loading", fnp

    cur.execute("""
DROP TABLE IF EXISTS {table};

CREATE TABLE {table}(
id serial PRIMARY KEY,
chrom varchar(31),
chromStart numeric,
chromEnd numeric,
gene text
);
""".format(table=table_name))

    ggff = Genes(args.file, "gtf")

    outF = StringIO.StringIO()
    for g in ggff.getGenes():
        outF.write("\t".join(g[:4]) + "\n")

    outF.seek(0)
    cur.copy_from(outF, table_name, ',',
                  columns=("chrom", "chromStart", "chromEnd", "gene"))

    cur.execute("""
CREATE INDEX {table}_idx01 ON {table}(gene);
""".format(table=table_name))

def setupAll(cur):
    setupAndCopy(cur, Dirs.GenomeFnp("gencode.m1/gencode.vM1.annotation.gtf.gz"),
                 "genes_mm9")
    #setupAndCopy(cur, Genome.mouse_gencode_m1_tss, "genes_hg19")
    #setupAndCopy(cur, Genome.mouse_gencode_m4_tss, "genes_mm10")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--local', action="store_true", default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    if args.local:
        dbs = DBS.localAnnotations()
    else:
        dbs = DBS.pgdsn("Annotations")
    dbs["application_name"] = os.path.realpath(__file__)

    import psycopg2.pool
    DBCONN = psycopg2.pool.ThreadedConnectionPool(1, 32, dbs)

    with getcursor(DBCONN, "main") as cur:
        if 0:
            setupAll(cur)

if __name__ == '__main__':
    main()
