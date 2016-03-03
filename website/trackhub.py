#!/usr/bin/env python

import os, sys, json
import StringIO

from helpers_trackhub import Track, PredictionTrack, BigGenePredTrack, BigWigTrack, officialVistaTrack, bigWigFilters

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from utils import Utils
from files_and_paths import Dirs
from web_epigenomes import WebEpigenome

BIB5 = "http://bib5.umassmed.edu/~purcarom/annotations_demo/"

eidToBigWigFileID = {"ENCSR172RHR" : "ENCFF301ECS",
                     "ENCSR179PIH" : "ENCFF138KVI",
                     "ENCSR196VDE" : "ENCFF803VYR",
                     "ENCSR292QBA" : "ENCFF329CUO",
                     "ENCSR312QVY" : "ENCFF592CFE",
                     "ENCSR337EDG" : "ENCFF849GKV",
                     "ENCSR358ESL" : "ENCFF989TMR",
                     "ENCSR367FCW" : "ENCFF147ZSW",
                     "ENCSR469VGZ" : "ENCFF142WSR",
                     "ENCSR488VEQ" : "ENCFF483YHV",
                     "ENCSR636NXY" : "ENCFF377GSQ",
                     "ENCSR661HDP" : "ENCFF128FJP",
                     "ENCSR666HFH" : "ENCFF479KSE",
                     "ENCSR687EAW" : "ENCFF091VKI",
                     "ENCSR742DUR" : "ENCFF385BSB",
                     "ENCSR793WUR" : "ENCFF957FAN",
                     "ENCSR935RRY" : "ENCFF877NOP",
                     "ENCSR959HKR" : "ENCFF617HTH",
                     "ENCSR064DGY" : "ENCFF478ATY",
                     "ENCSR080GQM" : "ENCFF337DHW",
                     "ENCSR157LYR" : "ENCFF096GTH",
                     "ENCSR196ENU" : "ENCFF079KZH",
                     "ENCSR293ORS" : "ENCFF749FRR",
                     "ENCSR335WME" : "ENCFF987TFU",
                     "ENCSR417TXZ" : "ENCFF937LBE",
                     "ENCSR425FLT" : "ENCFF588JVO",
                     "ENCSR464MQU" : "ENCFF686SHG",
                     "ENCSR581FAT" : "ENCFF247FWC",
                     "ENCSR641EME" : "ENCFF577TCZ",
                     "ENCSR765SJF" : "ENCFF884JJT",
                     "ENCSR831YAX" : "ENCFF354SIU",
                     "ENCSR871CGP" : "ENCFF167FIF",
                     "ENCSR929GXP" : "ENCFF350NWQ",
                     "ENCSR953LFI" : "ENCFF990QBR",
                     "ENCSR871KVM" : "ENCFF309IVZ",
                     "ENCSR446MUM" : "ENCFF077KNV",
                     "ENCSR319PWR" : "ENCFF722AND",
                     "ENCSR855ASN" : "ENCFF667HTC",
                     "ENCSR998KYQ" : "ENCFF986EOX",
                     "ENCSR767AJS" : "ENCFF061MOE",
                     "ENCSR791AJY" : "ENCFF788SQF",
                     "ENCSR723IXU" : "ENCFF174VMZ"
                     }

class TempWrap:
    def __init__(self, expID, fileID):
        self.expID = expID
        self.fileID = fileID
        self.url = os.path.join(BIB5, "data", expID, fileID + ".bigWig")
        self.file_status = "not known"

    def isBigWig(self):
        return self.url.endswith(".bigWig")

class TrackHub:
    def __init__(self, args, epigenomes, row):
        self.args = args
        self.epigenomes = epigenomes
        self.assembly = row[0]
        self.assays = row[1]
        self.tissue_ids = json.loads(row[2])
        self.loci = row[3]
        self.hubNum = row[4]

        self.priority = 1

    def Custom(self):
        lines = []
        #lines += ["browser hide all"]
        #lines += ["browser pack knownGene refGene ensGene"]
        #lines += ["browser dense snp128"]

        f = StringIO.StringIO()
        map(lambda line: f.write(line + "\n"), lines)

        return f.getvalue()

    def ParsePath(self, path):
        if not path:
            raise Exception("no path")

        if 1 == len(path):
            if path[0].startswith("hub_") and path[0].endswith(".txt"):
                return self.makeHub()
            if path[0].startswith("genomes_") and path[0].endswith(".txt"):
                return self.makeGenomes()
            return "ERROR"

        if 2 != len(path):
            raise Exception("path too long")

        if path[0] in ["hg19", "hg38", "mm9", "mm10"]:
            if path[0] == self.assembly:
                if path[1].startswith("trackDb_") and path[1].endswith(".txt"):
                    return self.makeTrackDb()

        raise Exception("invalid path")

    def makeHub(self):
        f = StringIO.StringIO()
        t = ""
        if self.args.debug:
            t += "debug "
        t += "ENCODE Encyclopedia Annotations " + self.assembly
        for r in [["hub", t],
                  ["shortLabel", t],
                  ["longLabel", t],
                  ["genomesFile", "genomes_{hubNum}.txt".format(hubNum=self.hubNum)],
                  ["email", "zhiping.weng@umassmed.edu"]]:
            f.write(" ".join(r) + "\n")
        return f.getvalue()

    def makeGenomes(self):
        return """genome\t{assembly}
trackDb\t{assembly}/trackDb_{hubNum}.txt""".format(assembly = self.assembly,
                                                   hubNum = self.hubNum)

    def makeTrackDb(self):
        epis = self.epigenomes.GetByAssemblyAndAssays(self.assembly, self.assays)
        epis = filter(lambda e: e.web_id() in self.tissue_ids, epis.epis)

        lines = []
        lines += [self.genes()]

        for wepi in sorted(epis, key=lambda e: e.epi.biosample_term_name):
            if "Both" == self.assays:
                lines += [self.predictionTrackHub(wepi)]
                lines += [self.compositeTrack(wepi)]
            for exp in wepi.exps():
                try:
                    lines += [self.trackhubExp(exp)]
                except:
                    if self.args.debug:
                        raise
                    pass

        if self.enableVistaTrack():
            lines += [self.vista()]
        lines += [self.phastcons()]

        lines = filter(lambda x: x, lines)

        f = StringIO.StringIO()
        map(lambda line: f.write(line + "\n"), lines)

        return f.getvalue()

    def enableVistaTrack(self):
        if "mm10" == self.assembly:
            for t in self.tissue_ids:
                if "11.5" in t:
                    return True
        return False

    def phastcons(self):
        if "mm9" == self.assembly:
            url = os.path.join(BIB5, "conservation", "mm9.phastCons30way.bw")
        if "mm10" == self.assembly:
            url =  "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw"
        if "hg19" == self.assembly:
            url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons100way/hg19.100way.phastCons.bw"

        desc = "phastCons"

        track = BigWigTrack(desc, self.priority, url, "0,255,0").track()
        self.priority += 1
        return track

    def genes(self):
        if "hg19" == self.assembly:
            return None

        byAssembly = {"mm9" : "Comprehensive M1",
                      "mm10" : "Comprehensive M8",
                      "hg19" : "Comprehensive 24"}
        desc = "GENCODE Genes " + byAssembly[self.assembly]

        byAssemblyURl = {"mm9" : os.path.join(BIB5, "genes", "gencode.vM1.annotation.bb"),
                         "mm10" : os.path.join(BIB5, "genes", "gencode.vM8.annotation.bb"),
                         "hg19" : os.path.join(BIB5, "genes", "gencode.v24.annotation.bb")}
        url = byAssemblyURl[self.assembly]

        track = BigGenePredTrack(desc, self.priority, url).track()
        self.priority += 1
        return track

    def vista(self):
        return officialVistaTrack(self.assembly)

    def predictionTrackHub(self, wepi):
        fnp = wepi.predictionFnp()
        if not os.path.exists(fnp):
            return None

        desc = Track.MakeDesc("candidate enhancers",
                              wepi.epi.age_display,
                              wepi.epi.biosample_term_name)

        url = os.path.join(BIB5, "Enhancer-Prediction-Tracks",
                           os.path.basename(fnp))

        track = PredictionTrack(desc, self.priority, url).track()
        self.priority += 1
        return track

    def trackhubExp(self, exp):
        url, name, color = self._getUrl(exp, False)

        desc = Track.MakeDesc(name, exp.age, exp.biosample_term_name)

        track = BigWigTrack(desc, self.priority, url, color).track()
        self.priority += 1
        return track

    def _getUrl(self, exp, norm):
        if not exp:
            return None, None, None

        bigWigs = bigWigFilters(self.assembly, exp.files)
        if not bigWigs:
            if "mm10" == self.assembly:
                bigWigs = [TempWrap(exp.encodeID,
                                    eidToBigWigFileID[exp.encodeID])]

        if not bigWigs:
            raise Exception("missing bigWigs for " + exp.encodeID)
        bigWig = bigWigs[0]

        url = bigWig.url
        if "mm10" == self.assembly and "released" != bigWig.file_status:
            url = os.path.join(BIB5, "data", bigWig.expID, bigWig.fileID + ".bigWig")

        if norm:
            if "mm9" == self.assembly or "mm10" == self.assembly:
                url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
            else:
                if bigWig.expID.startswith("EN"):
                    url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
                else:
                    url = os.path.join(BIB5, "roadmap_norm/consolidated/",
                                       bigWig.expID,
                                       bigWig.fileID + ".norm.bigWig")

        if exp.isH3K27ac():
            name = "H3K27ac Signal"
            color = "18,98,235"
        elif exp.isDNaseSeq():
            name = "DNase Signal"
            color = "255,121,3"
        else:
            raise Exception("unexpected exp")

        return url, name, color

    def compositeTrack(self, wepi):
        dnaseExp, h3k27acExp = wepi.exps()
        h3k27acUrl, h3k27acName, h3k27acColor = self._getUrl(h3k27acExp, True)
        dnaseUrl, dnaseName, dnaseColor = self._getUrl(dnaseExp, True)

        desc = wepi.web_title()
        descShort = "test"

        track = """
track composite{priority}
container multiWig
aggregate transparentOverlay
showSubtrackColorOnUi on
type bigWig 0 10.0
shortLabel {descShort}
longLabel {desc}
visibility full
priority {priority}
html examplePage

                track composite{priority}H3K27ac
                bigDataUrl {h3k27acUrl}
                shortLabel H3K27ac
                longLabel H3K27ac
                parent composite{priority}
                type bigWig
                color {h3k27acColor}

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
           h3k27acUrl = h3k27acUrl,
           h3k27acColor = h3k27acColor,
           dnaseUrl = dnaseUrl,
           dnaseColor = dnaseColor)

        self.priority += 1
        return track

    def showMissing(self, urlCache):
        wepis = self.epigenomes.GetByAssemblyAndAssays(self.assembly, self.assays)

        def checkUrl(url):
            if not url in urlCache:
                urlCache[url] = Utils.checkIfUrlExists(url)
            if urlCache[url]:
                return ""
            return url

        def checkExp(exp):
            u, _, _ = self._getUrl(exp, False)
            u = checkUrl(u)
            un, _, _ = self._getUrl(exp, True)
            un = checkUrl(un)
            return u, un

        for wepi in wepis.epis:
            dnaseExp = None
            h3k27acExp = None
            exps = wepi.exps()
            if "Both" == self.assays:
                dnaseExp, h3k27acExp = exps
            if "H3K27ac" == self.assays:
                h3k27acExp = exps[0]
            if "DNase" == self.assays:
                dnaseExp = exps[0]

            desc = wepi.web_title()
            dnaseUrl, dnaseUrlNorm = checkExp(dnaseExp)
            h3k27acUrl, h3k27acUrlNorm = checkExp(h3k27acExp)
            yield(desc, dnaseUrl, dnaseUrlNorm,
                  h3k27acUrl, h3k27acUrlNorm)

