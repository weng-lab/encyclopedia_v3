#!/usr/bin/env python

import os, sys, json, cherrypy, jinja2, argparse
import numpy as np
import uuid
import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.trackhub import TrackHub
from common.trackhub_washu import TrackHubWashu
from common.genes import LookupGenes
from common.dbsnps import dbSnps
from common.tables import DbTables
from common.session import Sessions
from common.db_trackhub import DbTrackhub
from common.db_bed_overlap import DbBedOverlap
from common.db_url_status import UrlStatusDB
from common.epigenome_stats import EpigenomeStats
from common.enums import AssayType
from common.ucsc_search import UcscSearch
from common.site_info import EnhancersSiteInfo

from models.enhancers.defaults import Defaults

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../metadata/utils'))
from utils import Utils
from templates import Templates

class TrackhubSite(object):
    def __init__(self, DBCONN, args, wepigenomes):
        self.args = args

        self.db = DbTrackhub(DBCONN)
        self.db_bed_overlap = DbBedOverlap(DBCONN)
        self.sessions = Sessions(DBCONN)
        self.dbSnps = dbSnps(DBCONN)
        self.genes = LookupGenes(DBCONN)
        self.wepigenomes = wepigenomes
        self.urlStatus = UrlStatusDB(DBCONN)

        viewDir = os.path.join(os.path.dirname(__file__), "../../views")
        self.templates = Templates(viewDir)

        self.host = "http://zlab-annotations.umassmed.edu/"
        if self.args.local:
            fnp = os.path.expanduser("~/.ws_host.txt")
            if os.path.exists(fnp):
                self.host = open(fnp).read().strip()

    @cherrypy.expose
    def index(self, *args, **params):
        return "index not defined"

    @cherrypy.expose
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def washu(self, *args, **params):
        uid = self.sessions.get(cherrypy.session.id)
        if not uid:
            uid = self.makeUid()
            cherrypy.session["uid"] = uid
            self.sessions.insert(cherrypy.session.id, uid)

        input_json = cherrypy.request.json

        us = UcscSearch(self.wepigenomes, self.db, self.dbSnps, self.genes,
                        self.host, self.args, input_json, uid)
        us.parse(self.siteInfo.site)
        url = us.configureWashuHubLink()

        if us.psb.userErrMsg:
            return { "err" : us.psb.userErrMsg }

        if self.args.debug:
            return {"inner-url" : url,
                    "html" : self.templates(self.siteInfo.site + "/ucsc",
                                            us = us,
                                            url = url)}
        return {"url" : url}

    @cherrypy.expose
    def trackhubCustom(self, *args, **params):
        cherrypy.response.headers['Content-Type'] = 'text/plain'

        uid = args[0]
        row = self.db.get(uid)
        if not row:
            raise Exception("uuid not found")

        th = TrackHub(self.args, self.wepigenomes, self.urlStatus, row)
        return th.Custom()

    @cherrypy.expose
    def trackhub(self, *args, **params):
        cherrypy.response.headers['Content-Type'] = 'text/plain'

        uid = args[0]
        row = self.db.get(uid)
        if not row:
            raise Exception("uuid not found")

        if AssayType.TargetGene == row["assayType"]:
            th = TrackHubTargetGene(self.args, self.wepigenomes, self.urlStatus, row)
        else:
            th = TrackHub(self.args, self.wepigenomes, self.urlStatus, row)

        path = args[1:]
        return th.ParsePath(path)

    @cherrypy.expose
    def trackhub_washu(self, *args, **params):
        cherrypy.response.headers['Content-Type'] = 'text/plain'

        uid = args[0]
        row = self.db.get(uid)
        if not row:
            raise Exception("uuid not found")

        th = TrackHubWashu(self.args, self.wepigenomes, self.urlStatus, row)

        path = args[1:]
        return th.ParsePath(path)

