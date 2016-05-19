#!/usr/bin/env python

from __future__ import print_function
import os, sys, json, psycopg2, argparse

from enums import AssayType

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from utils import Utils
from db_utils import getcursor, Databases
from tables import DbTables

class Sessions:
    def __init__(self, DBCONN, assay_type):
        self.DBCONN = DBCONN
        self.assay_type = assay_type

    @staticmethod
    def setupDB(DBCONN):
        with getcursor(DBCONN, "get") as curs:
            curs.execute("""
DROP TABLE IF EXISTS {table};
CREATE TABLE {table}
(id serial PRIMARY KEY,
uid text,
session_id text,
assay_type integer,
hub_num integer
) """.format(table = DbTables.sessions))
            print("recreated", DbTables.sessions, "table")

    def insert(self, session_id, uid, hub_num):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
INSERT INTO {table}
(session_id, uid, assay_type)
VALUES (
%(session_id)s,
%(uid)s,
%(assay_type)s,
%(hun_num)s
)""".format(table = DbTables.sessions), {"session_id" : session_id,
                                         "uid" : uid,
                                         "assay_type" : self.assay_type,
                                         "hub_num" : hub_num
})

    def get(self, session_id, hub_num):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
SELECT uid
FROM {table}
WHERE session_id = %(session_id)s
AND assay_type = %(assay_type)s
AND hub_num = %(hub_num)s
""".format(table = DbTables.sessions), {"session_id" : session_id,
                                        "assay_type" : self.assay_type,
                                        "hub_num" : hub_num
})
            uid = curs.fetchone()
            if uid:
                return uid[0]
            return None

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--local', action="store_true", default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    DBCONN = Databases.getAnnotationsDBCONN(args)

    Sessions.setupDB(DBCONN)

if __name__ == '__main__':
    main()
