#!/usr/bin/env python

import os, sys, json

from common.urls import BIB5

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from epigenome import Epigenomes
from exp_file import ExpFile
from files_and_paths import Dirs

UrlBase = os.path.join(BIB5, "Target-Gene-Prediction-Tracks/convert/")

class InteractingGeneExp:
    def __init__(self):
        self.encodeID = "InteractingGeneID"
        self.eid = self.encodeID
        self.assay_term_name = "InteractingGene"
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

class InteractingGeneEpigenome:
    def __init__(self, interacting_gene):
        self.assembly = "hg19"
        biosample_term_name = "GM12878"
        self.biosample_term_name = biosample_term_name
        self.biosample_term_id = ""
        self.organ_slims = ["blood"]
        self.tissue = "blood"
        self.biosample_type = "immortalized cell line"
        self.age_display = None
        self.eid = biosample_term_name

        self.interacting_gene = interacting_gene

    def hasInteractingGene(self):
        return True

    def InteractingGene(self):
        return [self.interacting_gene]

    def tadFnp(self, assays):
        path = Dirs.interactingGeneTracks
        if "TAD" == assays:
            fn = "{eid}_H3K4me3_predictions.bigBed".format(eid = self.eid)
        return os.path.join(path, fn)

class InteractingGeneMetadata:
    def __init__(self):
        self.epigenomes = Epigenomes("InteractingGene", "hg19")

        interacting_gene = InteractingGeneExp()
        epi = InteractingGeneEpigenome(interacting_gene)
        self.epigenomes.addEpigenome(epi)

        print "found", len(self.epigenomes), "epigenomes for Dekker hg19"

def main():
    r = InteractingGeneMetadata().epigenomes

if __name__ == '__main__':
    main()
