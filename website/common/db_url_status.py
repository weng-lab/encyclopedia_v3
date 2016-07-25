#!/usr/bin/env python

import os, sys, json, psycopg2, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from utils import Utils
from dbs import DBS
from db_utils import getcursor
from tables import DbTables

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

    adb = UrlStatusDB(DBCONN)
    adb.setupDB()

if __name__ == '__main__':
    main()
