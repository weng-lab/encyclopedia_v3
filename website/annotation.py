#!/usr/bin/env python

import os, sys, json, cherrypy, jinja2, argparse
import numpy as np
import uuid
import StringIO

from db import AnnotationDB
from session import Sessions
from ucsc_search import UcscSearch

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from utils import Utils
from trackhub import TrackHub
from web_epigenomes import WebEpigenomesLoader
from dbsnps import dbSnps
from defaults import Defaults
from parse_search_box import ParseSearchBox
from genes import LookupGenes

class Templates:
    def __init__(self, viewDir):
        self.views = jinja2.Environment(loader=jinja2.FileSystemLoader(viewDir))

    def __call__(self, t, **kwargs):
        if "title" not in kwargs:
            kwargs["title"] = ""
        if "meta" not in kwargs:
            kwargs["meta"] = []
        return self.views.get_template(t+".html").render(kwargs)

class AnnotationSearchUcsc(object):
    def __init__(self, DBCONN, args):
        self.args = args

        self.db = AnnotationDB(DBCONN)
        self.sessions = Sessions(DBCONN)
        self.dbSnps = dbSnps(DBCONN)
        self.genes = LookupGenes(DBCONN)

        viewDir = os.path.join(os.path.dirname(__file__), "views")
        self.templates = Templates(viewDir)

        self.epigenomes = WebEpigenomesLoader(self.args)

        self.defaults = Defaults()

        self.urlCache = {}

        self.host = "http://zlab-annotations.umassmed.edu/"
        if self.args.local:
            fnp = os.path.expanduser("~/.ws_host.txt")
            if os.path.exists(fnp):
                self.host = open(fnp).read().strip()

    @cherrypy.expose
    def default(self, *args, **params):
        return self.templates("index",
                              epigenomes = self.epigenomes,
                              defaults = self.defaults)

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

        us = UcscSearch(self.epigenomes, self.db, self.dbSnps, self.genes,
                        self.host, self.args, input_json, uid)
        us.parse()
        url = us.configureUcscHubLink()

        if us.psb.userErrMsg:
            return { "err" : us.psb.userErrMsg }

        if self.args.debug:
            return {"html" : self.templates("ucsc",
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

        th = TrackHub(self.args, self.epigenomes, row)
        return th.Custom()

    @cherrypy.expose
    def trackhub(self, *args, **params):
        cherrypy.response.headers['Content-Type'] = 'text/plain'

        uid = args[0]
        row = self.db.get(uid)
        if not row:
            raise Exception("uuid not found")

        th = TrackHub(self.args, self.epigenomes, row)

        path = args[1:]
        return th.ParsePath(path)

    @cherrypy.expose
    @cherrypy.tools.json_in()
    @cherrypy.tools.json_out()
    def bedsInRange(self, *args, **params):
        input_json = cherrypy.request.json
        ret = None

        try:
            psb = ParseSearchBox(self.epigenomes, self.dbSnps, self.genes, input_json)
            coord = psb.search()
            if coord:
                expIDs = self.db.findBedOverlap(psb.assembly,
                                                coord.chrom,
                                                coord.start,
                                                coord.end)
                ret = self.epigenomes.getWebIDsFromExpIDs(psb.assembly,
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
            return self.templates("missing_list")
        row = [args[0], args[1], "[]", "loci", "hubNum"]
        th = TrackHub(self.args, self.epigenomes, row)
        missing = th.showMissing(self.urlCache)
        return self.templates("missing",
                              missing = missing)
