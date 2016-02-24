#!/usr/bin/env python

import os, sys, cherrypy, json, argparse

import psycopg2, psycopg2.pool

from annotation.annotation import AnnotationSearchUcsc

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from dbs import DBS

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--local', action="store_true", default=False)
    parser.add_argument('--debug', action="store_true", default=False)
    parser.add_argument('--dev', action="store_true", default=True)
    parser.add_argument('--port', default=9191)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    staticDir = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                              "annotation/views/static"))
    cacheDir = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                             "annotation/cache"))
    asu_config = {
        '/': {
            'tools.sessions.on' : True,
            'tools.sessions.timeout' : 60000,
            'tools.sessions.storage_type' : "file",
            'tools.sessions.storage_path' : cacheDir
            },
        '/static': {
            'tools.staticdir.on': True,
            'tools.staticdir.dir': staticDir
            }
        }

    if args.local:
        dbs = DBS.localAnnotations()
    else:
        dbs = DBS.pgdsn("Annotations")
    dbs["application_name"] = os.path.realpath(__file__)
    DBCONN = psycopg2.pool.ThreadedConnectionPool(1, 32, **dbs)

    cherrypy.tree.mount(AnnotationSearchUcsc(DBCONN, args), '/',
                        config=asu_config)

    if args.dev:
        cherrypy.config.update({'server.environment': "development", })
    cherrypy.config.update({'server.socket_host': '0.0.0.0', })
    cherrypy.config.update({'server.socket_port': int(args.port), })

    if not args.local:
        # fend off harassment by the cluster
        cherrypy.config.update({'server.socket_queue_size': 512})
        cherrypy.config.update({'server.thread_pool': 30})

    cherrypy.engine.start()
    cherrypy.engine.block()

if __name__ == '__main__':
    sys.exit(main())
