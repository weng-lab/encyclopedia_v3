#!/usr/bin/env python

import os, sys, json
import StringIO

from urls import BIB5

from helpers_trackhub import Track, PredictionTrack, BigGenePredTrack, BigWigTrack, officialVistaTrack, bigWigFilters, TempWrap
from enums import AssayType

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.colors_trackhub import PredictionTrackhubColors, EncodeTrackhubColors, OtherTrackhubColors
from common.site_info import SiteInfos

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from utils import Utils
from files_and_paths import Dirs

class TrackHubWashu:
    def __init__(self, args, epigenomes, urlStatus, row):
        self.assembly = row["assembly"]
        self.assays = row["assays"]
        self.tissue_ids = json.loads(row["tissues"])
        self.loci = row["loci"]
        self.assayType = row["assayType"]
        self.hubNum = row["hubNum"]

        self.args = args
        self.epigenomes = epigenomes[self.assayType]
        self.urlStatus = urlStatus
        self.histMark = SiteInfos[self.assayType].histMark

        self.priority = 1

    def ParsePath(self, path):
        if not path:
            raise Exception("no path")

        if 2 != len(path):
            raise Exception("invalid path length")

        if path[0] in ["hg19", "hg38", "mm10"]:
            if path[0] == self.assembly:
                if path[1].startswith("trackDb_") and path[1].endswith(".json"):
                    return self.makeTrackDb()

        raise Exception("invalid path")

    def makeTrackDb(self):
        epis = self.epigenomes.GetByAssemblyAndAssays(self.assembly, self.assays)
        epis = filter(lambda e: e.web_id() in self.tissue_ids, epis.epis)

        ret = []
        #ret += [self.genes()]

        for wepi in sorted(epis, key=lambda e: e.epi.biosample_term_name):
            if self.assays.startswith("BothDNaseAnd"):
                ret += [self.predictionTrackHub(wepi)]
                #ret += [self.compositeTrack(wepi)]
            for exp in wepi.exps():
                try:
                    ret += [self.trackhubExp(exp)]
                except:
                    if self.args.debug:
                        raise
                    pass

#        if self.enableVistaTrack():
#            ret += [self.vista()]
        ret += [self.phastcons()]

        ret = filter(lambda x: x, ret)

        f = StringIO.StringIO()
        return json.dumps(ret, f)

    def enableVistaTrack(self):
        if "mm10" == self.assembly:
            for t in self.tissue_ids:
                if "11.5" in t:
                    return True
        return False

    def phastcons(self):
        if "mm10" == self.assembly:
            url =  "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw"
        if "hg19" == self.assembly:
            url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw"

        desc = "phastCons"

        color = OtherTrackhubColors.Conservation.rgb
        track = BigWigTrack(desc, self.priority, url, color).track_washu()
        self.priority += 1
        return track

    def genes(self):
        if "hg19" == self.assembly:
            return None

        byAssembly = {"mm10" : "Comprehensive M8",
                      "hg19" : "Comprehensive 24"}
        desc = "GENCODE Genes " + byAssembly[self.assembly]

        byAssemblyURl = {"mm10" : os.path.join(BIB5, "genes", "gencode.vM8.annotation.bb"),
                         "hg19" : os.path.join(BIB5, "genes", "gencode.v24.annotation.bb")}
        url = byAssemblyURl[self.assembly]

        track = BigGenePredTrack(desc, self.priority, url).track_washu()
        self.priority += 1
        return track

    def vista(self):
        return officialVistaTrack(self.assembly)

    def predictionTrackHub(self, wepi):
        fnp = wepi.predictionFnp()
        if not os.path.exists(fnp):
            return None

        if AssayType.Enhancer == self.assayType:
            descBase = "enhancer-like"
            url = os.path.join(BIB5,
                               Dirs.enhancerTracksBase,
                               "washu",
                               os.path.basename(fnp).replace(".bigBed",
                                                             ".bed.gz"))
        elif AssayType.Promoter == self.assayType:
            descBase = "promoter-like"
            url = os.path.join(BIB5,
                               Dirs.promoterTracksBase,
                               "washu",
                               os.path.basename(fnp).replace(".bigBed",
                                                             ".bed.gz"))

        desc = Track.MakeDesc(descBase,
                              wepi.epi.age_display,
                              wepi.epi.biosample_term_name)

        track = PredictionTrack(desc, self.priority, url).track_washu()
        self.priority += 1
        return track

    def trackhubExp(self, exp):
        url, name, color = self._getUrl(exp, False)

        desc = Track.MakeDesc(name, exp.age, exp.biosample_term_name)

        track = BigWigTrack(desc, self.priority, url, color).track_washu()
        self.priority += 1
        return track

    def _getUrl(self, exp, norm):
        if not exp:
            return None, None, None

        assay = "DNase"
        if exp.isH3K27ac():
            assay = "H3K27ac"
        elif exp.isH3K4me3():
            assay = "H3K4me3"

        bigWigs = bigWigFilters(self.assembly, exp.files)

        if not bigWigs:
            raise Exception("missing bigWigs for " + exp.encodeID)
        bigWig = bigWigs[0]

        url = bigWig.url
        if 1: #self.urlStatus.find(url) and not self.urlStatus.get(url):
            url = os.path.join(BIB5, "data", bigWig.expID,
                               bigWig.fileID + ".bigWig")

        if norm:
            if "mm10" == self.assembly:
                url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
            else:
                if bigWig.expID.startswith("EN"):
                    url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
                else:
                    url = os.path.join(BIB5, "roadmap_norm/consolidated/",
                                       bigWig.expID,
                                       bigWig.fileID + '-' + assay + ".fc.signal.norm.bigWig")

        if exp.isH3K27ac():
            name = "H3K27ac Signal"
            color = EncodeTrackhubColors.H3K27ac_Signal.rgb
        elif exp.isH3K4me3():
            name = "H3K4me3 Signal"
            color = EncodeTrackhubColors.H3K4me3_Signal.rgb
        elif exp.isDNaseSeq():
            name = "DNase Signal"
            color = EncodeTrackhubColors.DNase_Signal.rgb
        else:
            raise Exception("unexpected exp")

        return url, name, color

    def compositeTrack(self, wepi):
        dnaseExp, histoneExp = wepi.exps()
        histoneUrl, histoneName, histoneColor = self._getUrl(histoneExp, True)
        dnaseUrl, dnaseName, dnaseColor = self._getUrl(dnaseExp, True)

        desc = wepi.web_title()
        descShort = desc

        track = """
track composite{priority}
container multiWig
aggregate transparentOverlay
showSubtrackColorOnUi on
type bigWig 0 50.0
maxHeightPixels 128:32:8
shortLabel {descShort}
longLabel {desc}
visibility full
priority {priority}
html examplePage

                track composite{priority}Histone
                bigDataUrl {histoneUrl}
                shortLabel {histone}
                longLabel {histone}
                parent composite{priority}
                type bigWig
                color {histoneColor}

                track composite{priority}DNase
                bigDataUrl {dnaseUrl}
                shortLabel DNase
                longLabel DNase
                parent composite{priority}
                type bigWig
                color {dnaseColor}
""".format(priority = self.priority,
           descShort = descShort,
           desc = desc,
           histoneUrl = histoneUrl,
           histoneColor = histoneColor,
           dnaseUrl = dnaseUrl,
           dnaseColor = dnaseColor,
           histone = self.histMark)

        self.priority += 1
        return track
