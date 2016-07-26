#!/usr/bin/env python

import os, sys

from enums import AssayType

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils'))
from epigenome import Epigenomes
from helpers_metadata import ExpFile
from files_and_paths import Dirs

UrlBase = "http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/"

class RoadmapExp:
    def __init__(self, eid, assay_term_name, biosample_term_name,
                 tissue, biosample_type, files, histMark):
        self.encodeID = eid
        self.eid = eid
        self.assay_term_name = assay_term_name
        self.biosample_term_name = biosample_term_name
        self.tissue = tissue
        self.biosample_type = biosample_type
        self.age = None
        self.histMark = histMark

        self.files = None
        if files.strip():
            fn = eid + '-' + self.histMark + ".fc.signal.bigwig"
            if "DNase-seq" == assay_term_name:
                fn = eid + "-DNase.fc.signal.bigwig"

            expF = ExpFile.fromRoadmap(eid, self.assay_term_name)
            expF.url = os.path.join(UrlBase, fn)
            expF.output_type = "fold change over control"
            expF.file_type = "bigWig"

            self.files = [expF]

    def isH3K4me3(self):
        return "H3K4me3" == self.assay_term_name

    def isH3K27ac(self):
        return "H3K27ac" == self.assay_term_name

    def isDNaseSeq(self):
        return "DNase-seq" == self.assay_term_name

    def getIDRnarrowPeak(self, args = None):
        fn = "{eid}-{histMark}.narrowPeak.gz".format(eid=self.eid, histMark=self.histMark)
        if self.isDNaseSeq():
            fn = "{eid}-DNase.hotspot.fdr0.01.peaks.bed.gz".format(eid=self.eid)
        bedFnp = os.path.join("/project/umw_zhiping_weng/0_metadata/roadmap/data/consolidated", self.eid, fn)
        return bedFnp, "hg19"

    def getSingleBigWigSingleFnp(self, args = None):
        fn = "{eid}-{histMark}.fc.signal.bigwig".format(eid=self.eid, histMark=self.histMark)
        if self.isDNaseSeq():
            fn = "{eid}-DNase.fc.signal.bigwig".format(eid=self.eid)
        bedFnp = os.path.join("/project/umw_zhiping_weng/0_metadata/roadmap/data/consolidated", self.eid, fn)
        return bedFnp, "hg19"

class RoadmapEpigenome:
    def __init__(self, eid, biosample_term_name, tissue, biosample_type,
                 DNase, histone, histMark, assayType):
        self.assembly = "hg19"
        self.biosample_term_name = biosample_term_name
        self.biosample_term_id = ""
        self.organ_slims = [tissue.lower()]
        self.tissue = tissue
        self.biosample_type = biosample_type
        self.age_display = None
        self.eid = eid
        self.histMark = histMark
        self.assayType = assayType

        self.DNaseExp = None
        if DNase.files:
            self.DNaseExp = DNase

        self.histoneExp = None
        if histone.files:
            self.histoneExp = histone

    def hasDNase(self):
        return self.DNaseExp

    def hasH3K27ac(self):
        if "H3K27ac" == self.histMark:
            return self.histoneExp
        raise Exception("no H3K27ac")

    def hasH3K4me3(self):
        if "H3K4me3" == self.histMark:
            return self.histoneExp
        raise Exception("no H3K4me3")

    def DNase(self):
        return filter(lambda x: x, [self.DNaseExp])

    def Histone(self, histMark):
        if "H3K27ac" == histMark:
            return self.H3K27ac()
        if "H3K4me3" == histMark:
            return self.H3K4me3()
        raise Exception("unknown histMark " + histMark)

    def H3K27ac(self):
        if "H3K27ac" == self.histMark:
            return filter(lambda x: x, [self.histoneExp])
        raise Exception("no H3K27ac")

    def hasBothDNaseAndH3K27ac(self):
        if "H3K27ac" == self.histMark:
            return self.hasDNase() and self.hasH3K27ac()
        raise Exception("no H3K27ac")

    def H3K4me3(self):
        if "H3K4me3" == self.histMark:
            return filter(lambda x: x, [self.histoneExp])
        raise Exception("no H3K4me3")

    def hasBothDNaseAndH3K4me3(self):
        if "H3K4me3" == self.histMark:
            return self.hasDNase() and self.hasH3K4me3()
        raise Exception("no H3K4me3")

    def predictionFnp(self):
        if self.assayType == AssayType.Enhancer:
            return self.enhancerLikeFnp()
        if self.assayType == AssayType.Promoter:
            return self.promoterLikeFnp()
        raise Exception("unknown assay type")

    def enhancerLikeFnp(self, assays, DNase, H3K27ac):
        path = Dirs.enhancerTracks
        if "H3K27ac" == assays:
            fn = "{eid}_H3K27ac_predictions.bigBed".format(eid = self.eid)
        elif "DNase" == assays:
            fn = "{eid}_DNase_predictions.bigBed".format(eid = self.eid)
        elif "BothDNaseAndH3K27ac" == assays:
            fn = "{eid}_predictions.bigBed".format(eid = self.eid)
        else:
            raise Exception("fell through")
        return os.path.join(path, fn)

    def promoterLikeFnp(self, assays, DNase, H3K27ac):
        path = Dirs.promoterTracks
        if "H3K4me3" == assays:
            fn = "{eid}_H3K4me3_predictions.bigBed".format(eid = self.eid)
        elif "DNase" == assays:
            fn = "{eid}_DNase_predictions.bigBed".format(eid = self.eid)
        elif "BothDNaseAndH3K4me3" == assays:
            fn = "{eid}_predictions.bigBed".format(eid = self.eid)
        else:
            raise Exception("fell through")
        return os.path.join(path, fn)

class RoadmapMetadata:
    def __init__(self, histMark, assayType):
        self.histMark = histMark

        fnp = os.path.realpath(os.path.join(os.path.dirname(__file__), "roadmap.tsv"))
        with open(fnp) as f:
            data = [line.rstrip().split('\t') for line in f]
        headers = data[:3]
        data = data[3:]

        if 0:
            print headers[0][1]
            print headers[0][14]
            print headers[1][33]
            print headers[1][35]

            print data[0][1] # EID
            print data[0][2] # order
            print data[0][14] # biosample_term_name approx....
            print data[0][33] # H3K27ac tag align files
            print data[0][35] # DNase tag align files

        self.epigenomes = Epigenomes("ROADMAP", "hg19")

        for r in data:
            order = int(r[2])
            if order > 111: # exclude ENCODE2
                break
            DNase = RoadmapExp(r[1], "DNase-seq", r[14], r[16], r[17], r[35], histMark)
            histone = RoadmapExp(r[1], self.histMark, r[14], r[16], r[17], r[33], histMark)
            epi = RoadmapEpigenome(r[1], r[14], r[16], r[17], DNase, histone, histMark, assayType)
            self.epigenomes.addEpigenome(epi)

        print "found", len(self.epigenomes), "epigenomes for ROADMAP hg19"

def main():
    r = RoadmapMetadata().epigenomes

if __name__ == '__main__':
    main()
