import sys, os, cherrypy

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils'))
from utils import Utils
from templates import Templates

class VersionFour:
    def __init__(self):
        viewDir = os.path.join(os.path.dirname(__file__), "../views")
        self.templates = Templates(viewDir)

    @cherrypy.expose
    def index(self, *args, **params):
        raise cherrypy.HTTPRedirect("http://zlab-annotations-v4.umassmed.edu")
