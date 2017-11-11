#!/usr/bin/env python

import os
import sys

from parse_search_box import ParseSearchBox

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../metadata/utils'))
from files_and_paths import Dirs


class UcscSearch:
    def __init__(self, epigenomes, db, dbsnps, genes, host, args, params, uid):
        self.epigenomes = epigenomes
        self.db = db
        self.dbsnps = dbsnps
        self.genes = genes
        self.host = host
        self.args = args
        self.params = params
        self.uid = uid
        self.coord = None

    def Coord(self):
        if self.coord:
            return str(self.coord)
        return "None"

    def parse(self, siteInfo):
        try:
            self.psb = ParseSearchBox(self.epigenomes, self.dbsnps, self.genes, self.params)
            self.coord = self.psb.search()
            self.hubNum = self.db.insertOrUpdate(siteInfo.assayType,
                                                 self.psb.assembly,
                                                 self.psb.assays,
                                                 self.psb.tissue_ids,
                                                 self.psb.loci,
                                                 self.uid)
        except:
            raise
            if self.args.debug:
                raise
            pass

    def ucscParams(self):
        if self.coord:
            ucscParams = ["db=" + self.psb.assembly,
                          "position=" + str(self.coord)]
        else:
            # snp or gene
            if self.psb.assembly in ["hg19", "hg38"]:
                org = "Human"
            elif self.psb.assembly in ["mm10"]:
                org = "Mouse"
            else:
                raise Exception("unknown assembly")
            ucscParams = ["clade=mammal",
                          "org=" + org,
                          "db=" + self.psb.assembly,
                          "position=" + self.psb.loci,
                          "hgt.positionInput=" + self.psb.loci,
                          "hgt.suggestTrack=knownGene",
                          "Submit=submit"]
        if 0:
            customUrl = os.path.join(self.host,
                                     "trackhub/trackhub",
                                     "trackhubCustom",
                                     self.uid,
                                     str(self.hubNum))
            ucscParams.append("hgt.customText=" + customUrl)
        if 0:
            ucscParams = ["udcTimeout=1"] + ucscParams
        return ucscParams

    def configureUcscHubLink(self):
        ucscParams = self.ucscParams()

        urlBase = "https://genome.ucsc.edu/cgi-bin/hgTracks?"

        self.trackhubUrl = os.path.join(self.host,
                                        "trackhub/trackhub",
                                        self.uid,
                                        "hub_{hubNum}.txt".format(hubNum=self.hubNum))
        ucscParams.append("hubClear=" + self.trackhubUrl)

        self.trackdbUrl = os.path.join(self.host,
                                       "trackhub/trackhub",
                                       self.uid,
                                       self.psb.assembly,
                                       "trackDb_{hubNum}.txt".format(hubNum=self.hubNum))

        url = urlBase + "&".join(ucscParams)
        return url

    def configureWashuHubLink(self):
        self.trackdbUrl = os.path.join(self.host,
                                       "trackhub/trackhub_washu",
                                       self.uid,
                                       self.psb.assembly,
                                       "trackDb_{hn}.json".format(hn=self.hubNum))

        urlBase = "http://epigenomegateway.wustl.edu/browser/"
        assembly = "?genome=" + self.psb.assembly
        trackhub = "&datahub=" + self.trackdbUrl
        coord = "&coordinate=" + str(self.coord)

        url = urlBase + assembly + trackhub + coord
        return url
