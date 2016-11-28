#!/usr/bin/env python

from __future__ import print_function

import os, sys, json, argparse
from collections import defaultdict
from itertools import groupby
import cPickle as pickle

import numpy as np

from natsort import natsorted

from roadmap import RoadmapMetadata
from ontology.ontology import Ontology
from enums import AssayType

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from metadataws import MetadataWS
from files_and_paths import Datasets
from epigenome import Epigenomes

class ColWrap:
    def __init__(self, pretty_age, ageDays, selectorName):
        self.pretty_age = pretty_age
        self.ageDays = ageDays
        self.selectorName = selectorName

    def __eq__(self, other):
        return (isinstance(other, self.__class__)
                and self.__dict__ == other.__dict__)

    def __ne__(self):
         return not self.__eq__(other)

    # allow correct set() operations
    # http://stackoverflow.com/a/10547584
    def __hash__(self):
        return hash((self.pretty_age, self.selectorName))

    def matchWepi(self, wepi):
        return wepi.pretty_age() == self.pretty_age and wepi.SelectorName() == self.selectorName

class WalkRow:
    def __init__(self, row):
        self.row = row

    def hasSingleEntry(self):
        return 1 == np.count_nonzero(self.row)

    def Walk(self):
        justOne = self.hasSingleEntry()

        for exps in self.row:
            if exps == 0:
                yield None, None, None
            else:
                for c in exps:
                    if justOne:
                        yield c.web_id(), c.web_title_single(), c
                    else:
                        yield c.web_id(), c.web_title(), c

class WebEpigenomesLoader:
    def __init__(self, args, siteInfo):
        self.args = args
        self.ontology = Ontology()
        self.histMark = siteInfo.histMark
        self.assayType = siteInfo.assayType
        self.byAssemblyAssays = defaultdict(lambda : defaultdict(None))

        allAssays = ["BothDNaseAnd" + self.histMark, self.histMark, "DNase"]

        # mouse
        m = MetadataWS(Datasets.all_mouse)
        epigenomes = m.chipseq_tf_annotations_mm10(date = "2016-05-01")
        for assays in allAssays:
            self._loadEpigenomes("mm10", epigenomes, assays)

        # human
        roadmap = RoadmapMetadata(self.histMark, self.assayType).epigenomes
        m = MetadataWS(Datasets.all_human)
        encodeHg19 = m.chipseq_tf_annotations_hg19(date = "2016-05-01")
        combined = Epigenomes("ROADMAP + ENCODE", "hg19")
        combined.epis = encodeHg19.epis + roadmap.epis
        for assays in allAssays:
            self._loadEpigenomes("hg19", combined, assays)

    def _loadEpigenomes(self, assembly, epigenomes, assays):
        epis = epigenomes.GetByAssays(assays)
        if not epis:
            return
        wepis = []
        for epi in epis:
            if AssayType.Promoter == self.assayType:
                if epi.isRoadmap:
                    if assays.startswith("Both"):
                        # no human ROADMAP dual-ranked predictions
                        continue
            we = WebEpigenome(self.args, epi, assays,
                              self.ontology, self.histMark, self.assayType)
            wepis.append(we)

        self.byAssemblyAssays[assembly][assays] = WebEpigenomes(
            self.args, assembly, assays, wepis)
        
    def GetByAssemblyAndAssays(self, assembly, assays):
        if assembly not in self.byAssemblyAssays:
            print("missing", assembly)
            print(self.byAssemblyAssays.keys())
            raise Exception("missing " + assembly)
        if assays not in self.byAssemblyAssays[assembly]:
            print("missing", assays)
            print(self.byAssemblyAssays[assembly].keys())
            raise Exception("missing " + assays + " for " + assembly)
        return self.byAssemblyAssays[assembly][assays]

    def Walk(self, assembly, assays):
        return self.GetByAssemblyAndAssays(assembly, assays).Walk()

    def Header(self, assembly, assays):
        return self.GetByAssemblyAndAssays(assembly, assays).Header()

    def SelectorNames(self):
        ret = []
        for assembly in ["mm10", "hg19"]:
            for assays in ["BothDNaseAnd" + self.histMark, self.histMark, "DNase"]:
                epis = self.GetByAssemblyAndAssays(assembly, assays)
                for epi in epis.epis:
                    ret.append(epi.SelectorName())
        return json.dumps(sorted(list(set(ret))))

    def getWebIDsFromExpIDs(self, assembly, expIDs):
        ret = {}
        total = 0
        for assays in ["BothDNaseAnd" + self.histMark, self.histMark, "DNase"]:
            epis = self.GetByAssemblyAndAssays(assembly, assays)
            ret[assays] = epis.getWebIDsFromExpIDs(expIDs)
            total += len(ret[assays])
        ret["total"] = total
        return ret

class WebEpigenome:
    def __init__(self, args, epi, assays, ontology, histMark, assayType):
        self.args = args
        self.epi = epi
        self.assays = assays
        self.ontology = ontology
        self.histMark = histMark
        self.assayType = assayType

        self.DNase = None
        self.histones = None

        if len(self.epi.DNase()) > 1:
            print(self.epi)
            for e in self.epi.DNase():
                print("\t", e)
            raise Exception("multiple DNase experiments found")
        if len(self.epi.Histone(self.histMark)) > 1:
            print(self.epi)
            for e in self.epi.Histone(self.histMark):
                print("\t", e)
            raise Exception("multiple " + self.histMark + " experiments found")

        if "BothDNaseAnd" + self.histMark == self.assays:
            self.DNase = epi.DNase()[0]
            self.histones = epi.Histone(self.histMark)[0]
        elif self.histMark == self.assays:
            self.histones = epi.Histone(self.histMark)[0]
        elif "DNase" == self.assays:
            self.DNase = epi.DNase()[0]
        else:
            raise Exception("unknown assay type " + self.assays)

    def web_id(self):
        if self.epi.age_display:
            s = "_".join([self.epi.biosample_term_name, self.epi.life_stage, self.epi.age_display])
        else:
            if "hg19" == self.epi.assembly:
                s = "_".join([self.epi.biosample_term_name, "select"])
            else:
                s = "_".join([self.epi.biosample_term_name, "other"])
        return s.lower().replace(' ', '_').replace('.', '_')

    def pretty_age(self):
        s = ""
        if self.epi.age_display:
            s = self.epi.life_stage + " " + self.epi.age_display
        if "13.5 week" == self.epi.age_display:
            s = "e13.5" # exception for d embryonic fibroblast 13.5 week
        if "embryonic" in s:
            s = s.replace("embryonic ", "e").replace(" day", "")
        elif "postnatal" in s:
            s = s.replace(" day", "").replace("postnatal ", "p")
        if "adult 8 week" == s:
            s = "a8w"
        if "adult 8-10 week" == s:
            s = "a8-10w"
        if "adult 24 week" == s:
            s = "a24w"
        if not s:
            if "hg19" == self.epi.assembly:
                return ""
            else:
                s = "other"
        return s

    def ageDays(self):
        pa = self.pretty_age()
        lookup = { "e10" : 10,
                   "e10.5" : 10.5,
                   "e11" : 11,
                   "e11.5" : 11.5,
                   "e12" : 12,
                   "e13.5" : 13.5,
                   "e14.5" : 14.5,
                   "e15.5" : 15.5,
                   "e16.5" : 16.5,
                   "e18.5" : 18.5,
                   "p0" : 22,
                   "p1" : 23,
                   "p7" : 29,
                   "a8w" : 22 + 8*7,
                   "a8-10w" : 22 + 9*7,
                   "a24w" : 22 + 24*7,
                   "other" : 1000 }
        if pa not in lookup:
            print("ERROR: missing ageDays", pa)
            return 999
        return lookup[pa]

    def SelectorName(self):
        if self.epi.age_display:
            s = self.epi.life_stage + "_" + self.epi.age_display
            return s.replace('.', '_').replace(" ", "_")
        if "hg19" == self.epi.assembly:
            return "select"
        return "other"

    def web_title(self):
        if self.epi.age_display:
            return self.pretty_age()
        return self.epi.biosample_term_name

    def web_title_single(self):
        return self.pretty_age()

    def isActive(self):
        if AssayType.Enhancer == self.assayType:
            return self.web_id() in ["midbrain_embryonic_11_5_day",
                                     "hindbrain_embryonic_11_5_day",
                                     "limb_embryonic_11_5_day",
                                     "neural_tube_embryonic_11_5_day",
                                     "cerebellum_adult_8_week",
                                     "primary_t_cells_from_peripheral_blood_select",
                                     "primary_natural_killer_cells_from_peripheral_blood_select",
                                     "fetal_thymus_select"]
        if AssayType.Promoter == self.assayType:
            return self.web_id() in ["midbrain_embryonic_11_5_day",
                                     "hindbrain_embryonic_11_5_day",
                                     "limb_embryonic_11_5_day",
                                     "neural_tube_embryonic_11_5_day",
                                     "cerebellum_adult_8_week",
                                     "gm12878_select"]
        raise Exception("unknown assayType")

    def exps(self):
        if self.histMark == self.assays:
            return [self.histones]
        if "DNase" == self.assays:
            return [self.DNase]
        if "BothDNaseAnd" + self.histMark == self.assays:
            return [self.DNase, self.histones]
        raise Exception("unknown assay type " + self.assays)

    def enhancerLikePredictionFnp(self):
        return self.epi.enhancerLikeFnp(self.assays, self.DNase, self.histones)

    def enhancerLikePredictionFnpExists(self):
        fnp = self.enhancerLikePredictionFnp()
        ret = os.path.exists(fnp)
        if self.args.debug and not ret:
            print("enhancerLikePredictionFnpExists", "missing", self.epi.assembly, os.path.basename(fnp))
            print("\t-->", fnp)
        return ret

    def promoterLikePredictionFnp(self):
        fnp = self.epi.promoterLikeFnp(self.assays, self.DNase, self.histones)
        return fnp

    def promoterLikePredictionFnpExists(self):
        fnp = self.promoterLikePredictionFnp()
        ret = os.path.exists(fnp)
        if self.args.debug and not ret:
            print("promoterLikePredictionFnpExists", "missing",
                  self.epi.assembly, os.path.basename(fnp))
            print("\t-->", fnp)
        return ret

    def predictionFnp(self):
        if self.assayType == AssayType.Enhancer:
            return self.enhancerLikePredictionFnp()
        if self.assayType == AssayType.Promoter:
            return self.promoterLikePredictionFnp()
        raise Exception("unknown assay type")

    def predictionFnpExists(self):
        if self.assayType == AssayType.Enhancer:
            return self.enhancerLikePredictionFnpExists()
        if self.assayType == AssayType.Promoter:
            return self.promoterLikePredictionFnpExists()
        raise Exception("unknown assay type")

    def webName(self):
        if self.args.debug:
            if self.predictionFnpExists():
                return "*"
        return ""

    def webExps(self):
        if self.args.debug:
            return self.exps()
        return []

    def tissue_type(self):
        return self.ontology.getTissue(self.epi).replace('_', ' ')

    def cell_type(self):
        return self.ontology.getCellType(self.epi)

class WebEpigenomes:
    def __init__(self, args, assembly, assays, epis):
        self.args = args
        self.assembly = assembly
        self.assays = assays
        self.epis = epis

        self.expIDtoWebId = {}
        for wepi in self.epis:
            for exp in wepi.exps():
                self.expIDtoWebId[exp.encodeID] = wepi.web_id()

        self.setupMatrix()

    def __len__(self):
        return len(self.epis)

    def getWebIDsFromExpIDs(self, expIDs):
        return [self.expIDtoWebId[k] for k in expIDs if k in self.expIDtoWebId]

    def Header(self):
        for idx, c in enumerate(self.header):
            yield idx, c.pretty_age, c.selectorName

    def Walk(self):
        m = natsorted(self.m, key = lambda x: (x[0].lower(), x[1].lower(), x[2]))
        for row in m:
            yield row[0], row[1], row[2], WalkRow(row[3:]).Walk()

    def adjust_biosample_term_name(self, b):
        if "embryonic facial prominence" == b:
            b = "embryonic facial<br>prominence"
        b  = b.replace("negative", "-").replace("positive", "+").replace("--", "-").replace("-+", "+")
        return b

    def setupMatrix(self):
        if "hg19" == self.assembly:
            return self.setupMatrixHuman()
        return self.setupMatrixMouse()

    def setupMatrixHuman(self):
        wepis = self.epis
        keyfunc = lambda x: (x.epi.biosample_term_name.lower(), x.pretty_age())
	wepis.sort(key=keyfunc)

        cols = [ColWrap("Tissue of origin", "", ""),
                ColWrap("Cell Type", "", ""),
                ColWrap("Biosample", "", ""),
                ColWrap("", "", "select")]
        self.m = []

        # header row
        self.header = cols

        for rowIdx, wepi in enumerate(wepis):
            rowTitle = wepi.epi.biosample_term_name
            if wepi.pretty_age() and not wepi.epi.isImmortalizedCellLine():
                rowTitle += " (" + wepi.pretty_age() + ')'
            row = [wepi.tissue_type(),
                   wepi.cell_type(),
                   rowTitle,
                   [wepi]]
            self.m.append(row)

    def setupMatrixMouse(self):
        rows = set()
        cols = set()

        epis = self.epis
        keyfunc = lambda x: x.epi.biosample_term_name
	epis.sort(key=keyfunc)

        for biosample_term_name, wepis in groupby(epis, keyfunc):
            rows.add(biosample_term_name)
            for wepi in wepis:
                cols.add(ColWrap(wepi.pretty_age(), wepi.ageDays(), wepi.SelectorName()))

        cols = natsorted(list(cols), key = lambda x: x.ageDays)
        cols = [ColWrap("Tissue of origin", "", ""),
                ColWrap("Cell Type", "", ""),
                ColWrap("Biosample", "", "")
                ] + cols
        for colIdx, cw in enumerate(cols):
            if cw.pretty_age == "other":
                break
        cols += [cols.pop(colIdx)] # push "other" column to end
        self.m = []

        # header row
        self.header = cols

        for biosample_term_name, wepis in groupby(epis, keyfunc):
            row = [0] * len(cols)
            for wepi in wepis:
                # could do http://stackoverflow.com/a/947184
                for colIdx, cw in enumerate(cols):
                    if cw.matchWepi(wepi):
                        break
                if 0 == row[colIdx]:
                    row[colIdx] = []
                row[colIdx].append(wepi)
            row[0] = wepi.tissue_type()
            row[1] = wepi.cell_type()
            row[2] = self.adjust_biosample_term_name(biosample_term_name)
            self.m.append(row)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    we = WebEpigenomesLoader(args)

if __name__ == '__main__':
    main()
