#!/usr/bin/env python

import os
import sys
import json
import cherrypy
import jinja2

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../metadata/utils'))
from utils import Utils
from templates import Templates


class TutorialSite(object):
    def __init__(self, DBCONN, args, staticDir):
        self.DBCONN = DBCONN
        self.args = args
        self.staticDir = staticDir

        viewDir = os.path.join(os.path.dirname(__file__), "../../views")
        self.templates = Templates(viewDir)

        self.host = "http://zlab-annotations.umassmed.edu/"
        if self.args.local:
            fnp = os.path.expanduser("~/.ws_host.txt")
            if os.path.exists(fnp):
                self.host = open(fnp).read().strip()

    @cherrypy.expose
    def index(self, *args, **params):
        return self.templates('tutorial' + "/index",
                              site="tutorial")
