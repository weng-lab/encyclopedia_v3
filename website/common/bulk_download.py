import os, sys, json, cherrypy, time
import uuid
import StringIO
import zipfile

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.parse_search_box import ParseSearchBox

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils'))
from utils import Utils

class BulkDownload:
    def __init__(self, sessions, wepigenomes, dbSnps, genes, host, siteInfo, staticDir):
        self.sessions = sessions
        self.wepigenomes = wepigenomes
        self.dbSnps =dbSnps
        self.genes = genes
        self.host = host
        self.siteInfo = siteInfo
        self.staticDir = staticDir
        
    def makeUid(self):
        return str(uuid.uuid4())

    def download(self):
        uid = self.sessions.get(cherrypy.session.id)
        if not uid:
            uid = self.makeUid()
            cherrypy.session["uid"] = uid
            self.sessions.insert(cherrypy.session.id, uid)

        input_json = cherrypy.request.json

        try:
            psb = ParseSearchBox(self.wepigenomes, self.dbSnps, self.genes, input_json)
            if psb.userErrMsg:
                return { "err" : us.psb.userErrMsg }
        except:
            raise
            return { "err" : "internal download error"}

        epis = self.wepigenomes.GetByAssemblyAndAssays(psb.assembly, psb.assays)
        epis = filter(lambda e: e.web_id() in psb.tissue_ids, epis.epis)

        timestr = time.strftime("%Y%m%d-%H%M%S")
        outFn = timestr + '-' + '-'.join([self.siteInfo.name, psb.assembly,
                                          psb.assays]) + ".v3.zip"
        outFnp = os.path.join(self.staticDir, "downloads", uid, outFn)
        Utils.ensureDir(outFnp)

        counter = 0
        with zipfile.ZipFile(outFnp, mode='w', compression = zipfile.ZIP_STORED) as a:
            for wepi in epis:
                fnp = wepi.predictionFnp().replace(".bigBed", ".bed.gz")
                if not os.path.exists(fnp):
                    print("missing", fnp)
                    continue
                a.write(fnp, arcname = os.path.basename(fnp))
                counter += 1
        print("wrote", outFnp, counter)

        url = os.path.join(self.host, "static", "downloads", uid, outFn)
        return {"url" : url}
