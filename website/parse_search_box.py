#!/usr/bin/env python

import os, sys

from coord import Coord

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from files_and_paths import Dirs

class ParseSearchBox:
    def __init__(self, epigenomes, dbsnps, params):
        self.epigenomes = epigenomes
        self.dbsnps = dbsnps
        self.params = params

        self.halfWindow = 7500

        self._parse()

    def _parse(self):
        self.assembly = self.params["assembly"]
        self.assays = self.params["assays"]

        key = self.assembly + self.assays
        self.tissue_ids = []
        if key in self.params:
            self.tissue_ids = sorted(self.params[key])

        self.loci = self.params["loci"].strip()

        return self.search()

    def search(self):
        coord = None
        if self.loci.lower().startswith("chr"):
            # coordinate
            coord = Coord.parse(self.loci)
        elif self.loci.lower().startswith("rs"):
            # SNP
            coord = self.parseSnp()
        elif self.loci.isdigit() and 1 == len(self.tissue_ids):
            # a ranked peak for a single selected tissue
            coord = self.getRankedPeakCoord()
        return coord

    def parseSnp(self):
        snps = self.dbsnps.lookup(self.assembly, self.loci)
        if not snps:
            return None

        if len(snps) > 1:
            # search on UCSC
            return None

        snp = snps[0]
        c = Coord(snp[0], snp[1], snp[2])
        c.resize(self.halfWindow)
        return c

    def getRankedPeakCoord(self):
        wepis = self.epigenomes.GetByAssemblyAndAssays(self.assembly, self.assays)
        wepis = filter(lambda e: e.web_id() in self.tissue_ids, wepis.epis)

        if 1 != len(wepis):
            raise Exception("too many tissues")

        wepi = wepis[0]

        fnp = wepi.predictionFnp().replace(".bigBed", ".bed")
        if not os.path.exists(fnp):
            raise Exception("file not found " + fnp)

        rank = int(self.loci) - 1

        if rank < 0:
            raise Exception("invalid (negative) rank")

        row = None
        with open(fnp) as f:
            # http://stackoverflow.com/a/2081880
            for i, line in enumerate(f):
                if i == rank:
                    row = line
                    break
                if i > rank:
                    break
        if not row:
            raise Exception("invalid (too large) rank")

        toks = row.rstrip().split('\t')
        c = Coord(toks[0], toks[1], toks[2])
        c.resize(self.halfWindow)
        return c

