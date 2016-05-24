#!/usr/bin/env python

import os, sys, json, psycopg2, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from utils import Utils
from dbs import DBS
from db_utils import getcursor
from tables import DbTables

class DbTrackhub:
    def __init__(self, DBCONN):
        self.DBCONN = DBCONN
        self.tableSearch = DbTables.search

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
uid text NOT NULL,
site text NOT NULL,
hubNum integer NOT NULL
        ) """.format(search = self.tableSearch))

    def insert(self, site, assembly, assays, tissues, loci, uid):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT MAX(hubNum) FROM {search}
WHERE uid = %(uid)s
AND site = %(site)s
""".format(search = self.tableSearch), {"uid" : uid,
                                        "site" : site})
            hubNum = curs.fetchone()[0]
            if not hubNum:
                hubNum = 0
            else:
                hubNum += 1

            curs.execute("""
INSERT INTO {search}
(assembly, assays, tissues, loci, uid, site, hubNum)
VALUES (
%(assembly)s,
%(assays)s,
%(tissues)s,
%(loci)s,
%(uid)s,
%(site)s,
%(hubNum)s
);
""".format(search = self.tableSearch), {"assembly" : assembly,
                                        "assays" : assays,
                                        "tissues" : json.dumps(tissues),
                                        "loci" : loci,
                                        "uid" : uid,
                                        "site" : site,
                                        "hubNum" : hubNum
                                        })
        return hubNum

    def get(self, uid, hubNum):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
SELECT assembly, assays, tissues, loci, hubNum, site
FROM {search}
WHERE uid = %(uid)s
AND hubNum = %(hubNum)s
""".format(search = self.tableSearch), {"uid" : uid,
                                        "hubNum" : hubNum})
            return curs.fetchone()

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

    adb = DbTrackhub(DBCONN)
    adb.setupDB()

if __name__ == '__main__':
    main()
