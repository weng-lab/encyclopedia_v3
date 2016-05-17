#!/usr/bin/env python

import os, sys, json, psycopg2, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from utils import Utils
from dbs import DBS
from db_utils import getcursor
from tables import DbTables

class AnnotationDB:
    def __init__(self, DBCONN, tableSearch):
        self.DBCONN = DBCONN
        self.tableSearch = tableSearch

    def setupDB(self):
        with getcursor(self.DBCONN, "setupDB") as curs:
            curs.execute("""
    DROP TABLE IF EXISTS {search};
    CREATE TABLE {search}
    (id serial PRIMARY KEY,
    assembly text,
    assays text,
    tissues text,
    loci text,
    uid text UNIQUE NOT NULL,
    hubNum integer
        ) """.format(search = self.tableSearch))

    def findBedOverlap(self, assembly, chrom, start, end):
        if assembly not in ["hg19", "mm10"]:
            return []
        with getcursor(self.DBCONN, "findBedOverlap") as curs:
            curs.execute("""
SELECT DISTINCT expID
FROM bedRanges{assembly}
WHERE chrom = %(chrom)s
AND startend && int4range(%(start)s, %(end)s)
""".format(assembly=assembly),
{"chrom" : chrom,
 "start" : start,
 "end" : end})
            return [x[0] for x in curs.fetchall()]

    def insertOrUpdate(self, assembly, assays, tissues, loci, uid):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT id FROM {search}
WHERE uid = %(uid)s
""".format(search = self.tableSearch), {"uid" : uid})
            if (curs.rowcount > 0):
                curs.execute("""
UPDATE {search}
SET
assembly = %(assembly)s,
assays = %(assays)s,
tissues = %(tissues)s,
loci = %(loci)s,
hubNum = hubNum + 1
WHERE uid = %(uid)s
RETURNING hubNum;
""".format(search = self.tableSearch), {"assembly" : assembly,
      "assays" : assays,
      "tissues" : json.dumps(tissues),
      "loci" : loci,
      "uid" : uid
})
                hubNum = curs.fetchone()[0]
            else:
                curs.execute("""
INSERT INTO {search}
(assembly, assays, tissues, loci, uid, hubNum)
VALUES (
%(assembly)s,
%(assays)s,
%(tissues)s,
%(loci)s,
%(uid)s,
%(hubNum)s
) RETURNING hubNum;
""".format(search = self.tableSearch), {"assembly" : assembly,
       "assays" : assays,
       "tissues" : json.dumps(tissues),
       "loci" : loci,
       "uid" : uid,
       "hubNum" : 0
})
                hubNum = curs.fetchone()[0]
        return hubNum

    def get(self, uid):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
SELECT assembly, assays, tissues, loci, hubNum
FROM {search}
WHERE uid = %(uid)s
""".format(search = self.tableSearch), {"uid" : uid})
            return curs.fetchone()

class UrlStatusDB:
    def __init__(self, DBCONN):
        self.DBCONN = DBCONN

    def setupDB(self):
        with getcursor(self.DBCONN, "setupDB") as curs:
            curs.execute("""
DROP TABLE IF EXISTS urlExists;
CREATE TABLE urlExists
(id serial PRIMARY KEY,
url text,
exists boolean
    ) """)

    def find(self, url):
        with getcursor(self.DBCONN, "findBedOverlap") as curs:
            curs.execute("""
SELECT exists
FROM urlExists
WHERE url = %(url)s
""", {"url" : url})
            return curs.rowcount > 0

    def get(self, url):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
SELECT exists
FROM urlExists
WHERE url = %(url)s
""", {"url" : url})
            if curs.rowcount > 0:
                return curs.fetchone()[0]
            return False

    def insertOrUpdate(self, url, exists):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT exists FROM urlExists
WHERE url = %(url)s
""", {"url" : url})
            if (curs.rowcount > 0):
                curs.execute("""
UPDATE urlExists
SET
exists = %(exists)s
WHERE url = %(url)s
RETURNING id;
""", {"url" : url,
      "exists" : exists
})
            else:
                curs.execute("""
INSERT INTO urlExists
(url, exists)
VALUES (
%(url)s,
%(exists)s
)""", {"url" : url, "exists" : exists})

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

    adb = AnnotationDB(DBCONN, DbTables.search_promoters)
    adb.setupDB()

if __name__ == '__main__':
    main()
