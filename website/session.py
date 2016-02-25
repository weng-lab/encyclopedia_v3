#!/usr/bin/env python

import os, sys, json, psycopg2, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils
from dbs import DBS
from db_utils import getcursor

def setupDB(cur):
    cur.execute("""
DROP TABLE IF EXISTS sessions;
CREATE TABLE sessions
(id serial PRIMARY KEY,
uid text,
session_id text
) """)

class Sessions:
    def __init__(self, DBCONN):
        self.DBCONN = DBCONN

    def insert(self, session_id, uid):
        with getcursor(self.DBCONN, "get") as curs:
            curs.execute("""
INSERT INTO sessions
(session_id, uid)
VALUES (
%(session_id)s,
%(uid)s
)""", {"session_id" : session_id,
       "uid" : uid
})

    def insertOrUpdate(self, session_id, uid):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT id FROM sessions
WHERE session_id = %(session_id)s
""", {"session_id" : session_id})
            if (curs.rowcount > 0):
                curs.execute("""
UPDATE sessions
SET
uid = %(uid)s
WHERE session_id = %(session_id)s
""", {"session_id" : session_id,
      "uid" : uid
})
            else:
                curs.execute("""
INSERT INTO sessions
(session_id, uid)
VALUES (
%(session_id)s,
%(uid)s
)""", {"session_id" : session_id,
       "uid" : uid
})

    def get(self, session_id):
        with getcursor(self.DBCONN, "insertOrUpdate") as curs:
            curs.execute("""
SELECT uid
FROM sessions
WHERE session_id = %(session_id)s
""", {"session_id" : session_id})
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


    if args.local:
        dbs = DBS.localAnnotations()
    else:
        dbs = DBS.pgdsn("Annotations")
    dbs["application_name"] = os.path.realpath(__file__)

    import psycopg2.pool
    DBCONN = psycopg2.pool.ThreadedConnectionPool(1, 32, dbs)

    with getcursor(DBCONN, "main") as cur:
        setupDB(cur)

if __name__ == '__main__':
    main()
