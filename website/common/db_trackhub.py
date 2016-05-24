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
assayType integer NOT NULL
) """.format(search = self.tableSearch))

    def get(self, uid):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
SELECT assembly, assays, tissues, loci, assayType
FROM {search}
WHERE uid = %(uid)s
""".format(search = self.tableSearch), {"uid" : uid})
            return curs.fetchone()

    def insertOrUpdate(self, assayType, assembly, assays, tissues, loci, uid):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT id FROM search
WHERE uid = %(uid)s
""", {"uid" : uid})
            if (curs.rowcount > 0):
                curs.execute("""
UPDATE search
SET
assembly = %(assembly)s,
assays = %(assays)s,
tissues = %(tissues)s,
loci = %(loci)s,
assayType = %(assayType)s
WHERE uid = %(uid)s
""", {"assembly" : assembly,
      "assays" : assays,
      "tissues" : json.dumps(tissues),
      "loci" : loci,
      "uid" : uid,
      "assayType" : assayType
})
            else:
                curs.execute("""
INSERT INTO search
(assembly, assays, tissues, loci, uid, assayType)
VALUES (
%(assembly)s,
%(assays)s,
%(tissues)s,
%(loci)s,
%(uid)s,
%(assayType)s
);
""", {"assembly" : assembly,
      "assays" : assays,
      "tissues" : json.dumps(tissues),
      "loci" : loci,
      "uid" : uid,
      "assayType" : assayType
})

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
