#!/usr/bin/env python

import os, sys, json

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from epigenome import Epigenomes
from helpers_metadata import ExpFile
from files_and_paths import Dirs

UrlBase = "http://bib5.umassmed.edu/~purcarom/hic-Dekker/tads/"

TADS = { "SK-N-DZ" : "ENCODE3-SKNDZ-HindIII__hg19__ucsc/ENCODE3-SKNDZ-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "SK-MEL-5" : "ENCODE3-SKMEL5-HindIII__hg19__ucsc/ENCODE3-SKMEL5-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "A549" : "ENCODE3-A549-HindIII__hg19__ucsc/ENCODE3-A549-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "LNCaP clone FGC" : "ENCODE3-LNCaP-HindIII__hg19__ucsc/ENCODE3-LNCaP-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "Panc1" : "ENCODE3-PANC1-HindIII__hg19__ucsc/ENCODE3-PANC1-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "NCI-H460" : "ENCODE3-NCIH460-HindIII__hg19__ucsc/ENCODE3-NCIH460-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "SK-N-MC" : "ENCODE3-SKNMC-HindIII__hg19__ucsc/ENCODE3-SKNMC-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "G401" : "ENCODE3-G401-HindIII__hg19__ucsc/ENCODE3-G401-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "Caki2" : "ENCODE3-CAKI2-HindIII__hg19__ucsc/ENCODE3-CAKI2-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "SJCRH30" : "ENCODE3-SJCRH30-HindIII__hg19__ucsc/ENCODE3-SJCRH30-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "T47D" : "ENCODE3-T470-HindIII__hg19__ucsc/ENCODE3-T470-HindIII__hg19__genome__C-40000-iced.tads.bed",
         "RPMI-7951" : "ENCODE3-RPMI7951-HindIII__hg19__ucsc/ENCODE3-RPMI7951-HindIII__hg19__genome__C-40000-iced.tads.bed"}

class DekkerExp:
    def __init__(self, biosample_term_name, info):
        self.encodeID = info["accession"]
        self.eid = self.encodeID
        self.assay_term_name = "HiC"
        self.biosample_term_name = biosample_term_name
        self.tissue = biosample_term_name
        self.biosample_type = biosample_term_name
        self.age = None

        expF = ExpFile.fromRoadmap(self.eid, self.assay_term_name)
        expF.url = os.path.join(UrlBase, TADS[biosample_term_name])
        expF.output_type = "TAD"
        expF.file_type = "bed"

        self.files = [expF]

    def getIDRnarrowPeak(self, args = None):
        return None, "hg19"

    def getSingleBigWigSingleFnp(self, args = None):
        return None, "hg19"

class DekkerEpigenome:
    def __init__(self, biosample_term_name, info, tad):
        self.assembly = "hg19"
        self.biosample_term_name = biosample_term_name
        self.biosample_term_id = ""
        self.organ_slims = biosample_term_name
        self.tissue = biosample_term_name
        self.biosample_type = biosample_term_name
        self.age_display = None
        self.eid = info["accession"]

        self.tad = tad

    def hasTAD(self):
        return True

    def TAD(self):
        return [self.tad]

    def tadFnp(self, assays):
        path = Dirs.promoterTracks
        if "TAD" == assays:
            fn = "{eid}_H3K4me3_predictions.bigBed".format(eid = self.eid)
        return os.path.join(path, fn)

class DekkerMetadata:
    def __init__(self):
        fnp = os.path.realpath(os.path.join(os.path.dirname(__file__), "hic.txt.json"))
        with open(fnp) as f:
            data = json.load(f)

        self.epigenomes = Epigenomes("Dekker", "hg19")

        for biosample_term_name, info in data.iteritems():
            accession = info["accession"]
            TAD = DekkerExp(biosample_term_name, info)
            epi = DekkerEpigenome(biosample_term_name, info, TAD)
            self.epigenomes.addEpigenome(epi)

        print "found", len(self.epigenomes), "epigenomes for Dekker hg19"

def main():
    r = DekkerMetadata().epigenomes

if __name__ == '__main__':
    main()
