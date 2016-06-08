#!/usr/bin/env python

import os, sys, json

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from epigenome import Epigenomes
from helpers_metadata import ExpFile
from files_and_paths import Dirs

UrlBase = "http://bib5.umassmed.edu/~purcarom/annotations_demo/Target-Gene-Prediction-Tracks/convert/"

class TargetGeneExp:
    def __init__(self):
        self.encodeID = "TargetGeneID"
        self.eid = self.encodeID
        self.assay_term_name = "TargetGene"
        self.biosample_term_name = "GM12878"
        self.tissue = "blood"
        self.biosample_type = "immortalized cell line"
        self.age = None

        self.files = []
        for a in ["eQTL", "HiC", "POL2-ChIAPet", "RAD21-ChIAPet"]:
            expF = ExpFile.fromRoadmap(self.eid, a)
            expF.url = os.path.join(UrlBase, a + ".gtf.genePred")
            expF.output_type = a
            expF.file_type = "genePred"
            self.files.append(expF)

    def getIDRnarrowPeak(self, args = None):
        return None, "hg19"

    def getSingleBigWigSingleFnp(self, args = None):
        return None, "hg19"

class TargetGeneEpigenome:
    def __init__(self, target_gene):
        self.assembly = "hg19"
        biosample_term_name = "GM12878"
        self.biosample_term_name = biosample_term_name
        self.biosample_term_id = ""
        self.organ_slims = ["blood"]
        self.tissue = "blood"
        self.biosample_type = "immortalized cell line"
        self.age_display = None
        self.eid = biosample_term_name

        self.target_gene = target_gene

    def hasTargetGene(self):
        return True

    def TargetGene(self):
        return [self.target_gene]

    def tadFnp(self, assays):
        path = Dirs.targetGeneTracks
        if "TAD" == assays:
            fn = "{eid}_H3K4me3_predictions.bigBed".format(eid = self.eid)
        return os.path.join(path, fn)

class TargetGeneMetadata:
    def __init__(self):
        self.epigenomes = Epigenomes("TargetGene", "hg19")

        target_gene = TargetGeneExp()
        epi = TargetGeneEpigenome(target_gene)
        self.epigenomes.addEpigenome(epi)

        print "found", len(self.epigenomes), "epigenomes for Dekker hg19"

def main():
    r = TargetGeneMetadata().epigenomes

if __name__ == '__main__':
    main()
