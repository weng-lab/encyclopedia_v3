#!/usr/bin/env python

import os, sys, json, cherrypy, jinja2, argparse
import numpy as np
import uuid
import StringIO

from ucsc_search import UcscSearch
from trackhub import TrackHub
from trackhub_washu import TrackHubWashu
from parse_search_box import ParseSearchBox

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.genes import LookupGenes
from common.dbsnps import dbSnps

from models.enhancers.session import Sessions
from models.enhancers.web_epigenomes import WebEpigenomesLoader
from models.enhancers.db import AnnotationDB, UrlStatusDB
from models.enhancers.defaults import Defaults
from models.enhancers.epigenome_stats import EpigenomeStats

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../metadata/utils'))
from utils import Utils
from templates import Templates

class EnhancersSite(object):
    def __init__(self, DBCONN, args):
        self.args = args

        self.site = "enhancers"
        self.db = AnnotationDB(DBCONN)
        self.sessions = Sessions(DBCONN)
        self.dbSnps = dbSnps(DBCONN)
        self.genes = LookupGenes(DBCONN)
        self.urlStatus = UrlStatusDB(DBCONN)
        self.wepigenomes = WebEpigenomesLoader(self.args)
        self.defaults = Defaults()
        self.epigenome_stats = EpigenomeStats(self.wepigenomes)

        viewDir = os.path.join(os.path.dirname(__file__), "../../views")
        self.templates = Templates(viewDir)

        self.host = "http://zlab-annotations.umassmed.edu/"
        if self.args.local:
            fnp = os.path.expanduser("~/.ws_host.txt")
            if os.path.exists(fnp):
                self.host = open(fnp).read().strip()
        self.host += self.site + "/"

    @cherrypy.expose
    def index(self, *args, **params):
        return self.templates(self.site + "/index",
                              epigenomes = self.wepigenomes,
                              defaults = self.defaults,
                              stats = self.epigenome_stats,
                              site = self.site)

    def makeUid(self):
        return str(uuid.uuid4())

    @cherrypy.expose
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def ucsc(self, *args, **params):
        uid = self.sessions.get(cherrypy.session.id)
        if not uid:
            uid = self.makeUid()
            cherrypy.session["uid"] = uid
            self.sessions.insertOrUpdate(cherrypy.session.id, uid)

        input_json = cherrypy.request.json

        us = UcscSearch(self.wepigenomes, self.db, self.dbSnps, self.genes,
                        self.host, self.args, input_json, uid)
        us.parse()
        url = us.configureUcscHubLink()

        if us.psb.userErrMsg:
            return { "err" : us.psb.userErrMsg }

        if self.args.debug:
            return {"inner-url" : url,
                    "html" : self.templates(self.site + "/ucsc",
                                            us = us,
                                            url = url)}
        return {"url" : url}

    @cherrypy.expose
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def washu(self, *args, **params):
        uid = self.sessions.get(cherrypy.session.id)
        if not uid:
            uid = self.makeUid()
            cherrypy.session["uid"] = uid
            self.sessions.insertOrUpdate(cherrypy.session.id, uid)

        input_json = cherrypy.request.json

        us = UcscSearch(self.wepigenomes, self.db, self.dbSnps, self.genes,
                        self.host, self.args, input_json, uid)
        us.parse()
        url = us.configureWashuHubLink()

        if us.psb.userErrMsg:
            return { "err" : us.psb.userErrMsg }

        if self.args.debug:
            return {"inner-url" : url,
                    "html" : self.templates(self.site + "/ucsc",
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

    @cherrypy.expose
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def bedsInRange(self, *args, **params):
        input_json = cherrypy.request.json
        ret = None

        try:
            psb = ParseSearchBox(self.wepigenomes, self.dbSnps, self.genes, input_json)
            coord = psb.search()
            if coord:
                expIDs = self.db.findBedOverlap(psb.assembly,
                                                coord.chrom,
                                                coord.start,
                                                coord.end)
                ret = self.wepigenomes.getWebIDsFromExpIDs(psb.assembly,
                                                          expIDs)
            if psb.userErrMsg:
                return { "err" : psb.userErrMsg }

            return {"ret" : ret}

        except:
            if self.args.debug:
                raise

        return {"err" : "Problem parsing coordinate"}

    @cherrypy.expose
    def missing(self, *args, **params):
        if not args:
            return self.templates(self.site + "/missing_list")
        row = [args[0], args[1], "[]", "loci", "hubNum"]
        th = TrackHub(self.args, self.wepigenomes, self.urlStatus, row)
        missing = th.showMissing()
        return self.templates(self.site + "/missing",
                              missing = missing)

    @cherrypy.expose
    def methods(self, *args, **params):
        return self.templates(self.site + "/methods",
                              stats = self.epigenome_stats,
                              site = self.site)
