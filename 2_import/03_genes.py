#!/usr/bin/env python

import os
import sys
import psycopg2
import argparse
import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from dbs import DBS
from files_and_paths import Dirs
from db_utils import getcursor
from get_tss import Genes


def setupAndCopy(cur, fnp, fileType, table_name):
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

    ggff = Genes(fnp, fileType)

    outF = StringIO.StringIO()
    for g in ggff.getGenes():
        outF.write("\t".join([str(x) for x in [g.chr_, g.start_, g.end_, g.genename_]]) + "\n")

    outF.seek(0)
    cur.copy_from(outF, table_name, '\t',
                  columns=("chrom", "chromStart", "chromEnd", "gene"))
    print "\t", fnp, cur.rowcount
    cur.execute("""
    CREATE INDEX {table}_idx01 ON {table}(lower(gene));
""".format(table=table_name))


def setupAll(cur):
    setupAndCopy(cur, Dirs.GenomeFnp("gencode.m4/gencode.vM4.annotation.gtf.gz"),
                 "gtf", "genes_mm10")
    setupAndCopy(cur, Dirs.GenomeFnp("gencode.v19/gencode.v19.annotation.gff3.gz"),
                 "gff", "genes_hg19")


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
    DBCONN = psycopg2.pool.ThreadedConnectionPool(1, 32, **dbs)

    with getcursor(DBCONN, "main") as cur:
        setupAll(cur)


if __name__ == '__main__':
    main()
