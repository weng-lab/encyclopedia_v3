#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse
import json
import hashlib
import gzip
from pyliftover import LiftOver

import submit
from submit import getvalidator, EncodeFile
from arguments import *
from validator import validate_epigenome_metadata

from helpers_submit import authenticateEncodeTxt, getFileValidator, submitAnnotationMetadata, submitFile

sys.path.append(os.path.join(os.path.dirname(__file__), '../utils'))
from files_and_paths import Dirs, Tools
from querydcc import QueryDCC
from utils import Utils
from cache_memcache import MemCacheWrapper


class SubmitHiCLiftOver:
    def __init__(self, args):
        self.args = args
        self.doLiftOver = LiftOver('hg19', 'hg38')

        self.lengths_orig = []
        self.lengths_filtered = []
        self.oldVsNew = []

    def splitStrCoordStr(self, raw):
        chrom = raw.split(':')[0]
        start = raw.split(':')[1].split('-')[0]
        end = raw.split(':')[1].split('-')[1]
        return "\t".join([chrom, start, end])

    def splitStrCoord(self, raw):
        chrom = raw.split(':')[0]
        start = raw.split(':')[1].split('-')[0]
        end = raw.split(':')[1].split('-')[1]
        return [chrom, int(start), int(end)]

    def wrapLiftover(self, debug, chrom, start, end, errMsg):
        lift_start = self.doLiftOver.convert_coordinate(chrom, start)
        if not lift_start:
            if debug:
                print(errMsg + " start", chrom, start)
            return None
        lift_start = lift_start[0]
        lift_end = self.doLiftOver.convert_coordinate(chrom, end)
        if not lift_end:
            if debug:
                print(errMsg + " end", chrom, end)
            return None
        lift_end = lift_end[0]
        if lift_start[0] != lift_end[0]:
            if debug:
                print(errMsg + " no longer same chrom", chrom, start, end, lift_start[0], lift_end[0])
            return None
        oldLen = end - start

        chromLift = lift_start[0]
        startLift = lift_start[1]
        endLift = lift_end[1]
        newLen = endLift - startLift

        if oldLen < 1:
            if debug:
                print(errMsg + " oldLen: negative!", chrom, start, end)
            return None
        if newLen < 1:
            if debug:
                print(errMsg + " newLen: negative!", chromLift, startLift, endLift)
            return None

        absDiff = abs(newLen - oldLen)

        return [chromLift, startLift, endLift, oldLen, newLen, absDiff]

    def coordToStr(self, c):
        return c[0] + ':' + str(c[1]) + '-' + str(c[2])

    def parseLine(self, line):
        # chr10   3240001 4120000 boundary.3|hg19|chr10:3240001-3280000___boundary.4|hg19|chr10:4080001-4120000   1.06090369391
        # [0chrom, 1start, 2end, 3mess, 4value]
        toks = line.split()
        leftCoord = toks[:3]
        leftCoord[1] = int(leftCoord[1])
        leftCoord[2] = int(leftCoord[2])
        mtoks = toks[3].split('|')
        midBoundaryLeft = mtoks[0]

        if 3 != len(mtoks):
            midBoundaryRight = mtoks[2].split('__')[1]

        midCoordRaw = mtoks[2].split('__')[0]
        midCoord = self.splitStrCoord(midCoordRaw)

        if 3 != len(mtoks):
            rightCoord = self.splitStrCoord(mtoks[-1])

        leftCoordLift = self.wrapLiftover(False, leftCoord[0], leftCoord[1], leftCoord[2], "left")
        if not leftCoordLift:
            return None
        self.lengths_orig.append([leftCoordLift[3], leftCoordLift[4]])

        if leftCoordLift[5] > 5000:
            if 0:
                print("skipping b/c of lengths change")
            return None

        midCoordLift = self.wrapLiftover(False, midCoord[0], midCoord[1], midCoord[2], "mid")
        if not midCoordLift:
            return None
        if midCoordLift[5] > 5000:
            return None
        if 3 != len(mtoks):
            rightCoordLift = self.wrapLiftover(False, rightCoord[0], rightCoord[1], rightCoord[2], "right")
            if not rightCoordLift:
                return None
            if rightCoordLift[5] > 5000:
                return None

        self.lengths_filtered.append([leftCoordLift[3], leftCoordLift[4]])

        if 3 != len(mtoks):
            mid = [midBoundaryLeft, "hg38-liftOver", self.coordToStr(midCoordLift) + '___' + midBoundaryRight,
                   "hg38-liftOver", self.coordToStr(rightCoordLift)]
        else:
            mid = [midBoundaryLeft, "hg38-liftOver", self.coordToStr(midCoordLift)]

        ret = "\t".join([str(x) for x in leftCoordLift[:3] + ['|'.join(mid)] + [toks[4]]])
        self.oldVsNew.append([line, ret])
        return ret

    def tmpFile(self, accession, assembly, prefix):
        return os.path.join("/home/mjp/tadsLiftOverHg19ToHg38",
                            assembly + "_liftOver_" + prefix + '_' + accession + ".bed.gz")

    def parseOutFile(self, accession, fnp):
        good = 0
        bad = 0
        with gzip.open(fnp) as f:
            with gzip.open(self.tmpFile(accession, 'hg38', 'point'), 'wb') as outF:
                for line in f:
                    newLine = self.parseLine(line)
                    if newLine:
                        outF.write(newLine + '\n')
                        good += 1
                    else:
                        bad += 1
        print("lifted:", accession, good, bad)

    def runLiftover(self):
        mc = MemCacheWrapper()
        qd = QueryDCC(cache=mc)
        url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=Hi-C&status=released"

        for exp in qd.getExps(url):
            for f in exp.getTADs():
                f.download()
                self.parseOutFile(f.fileID, f.fnp())

        fnp = "/home/mjp/tadsLiftOverHg19ToHg38/lengths_orig.tsv"
        with open(fnp, 'w') as f:
            for r in self.lengths_orig:
                f.write('\t'.join([str(x) for x in r]) + '\n')
        print("wrote", fnp)

        fnp = "/home/mjp/tadsLiftOverHg19ToHg38/lengths_filtered.tsv"
        with open(fnp, 'w') as f:
            for r in self.lengths_filtered:
                f.write('\t'.join([str(x) for x in r]) + '\n')
        print("wrote", fnp)

        fnp = "/home/mjp/tadsLiftOverHg19ToHg38/oldVsNew.tsv"
        with open(fnp, 'w') as f:
            for r in self.oldVsNew:
                f.write(r[0])
                f.write(r[1] + '\n')
        print("wrote", fnp)

    def fileJson(self, exp, f, fnp):
        return {
            "dataset": exp.encodeID,
            "file_format": "bed",
            "file_format_type": "bed3+",
            "file_size": os.path.getsize(fnp),
            "md5sum": Utils.md5(fnp),
            "output_type": f.output_type,
            "assembly": "GRCh38",
            "award": "/awards/U41HG007000/",
            "lab": "/labs/zhiping-weng/",
            "derived_from": [f.fileID],
            "submitted_file_name": fnp,
            "aliases": ["zhiping-weng:hic-tad-hg38-liftOver-" + f.fileID]
        }

    def submitFile(self, exp, f):
        fileAccession = f.fileID
        fnp = self.tmpFile(fileAccession, 'hg38', 'point')
        j = self.fileJson(exp, f, fnp)
        print(j)
        submitFile(self.args, j)

    def runSubmit(self):
        authenticateEncodeTxt(self.args)

        mc = MemCacheWrapper()
        qd = QueryDCC(cache=mc)
        url = "https://www.encodeproject.org/search/?type=Experiment&assay_title=Hi-C&status=released"

        for exp in qd.getExps(url):
            for f in exp.getTADs():
                f.download()
                self.submitFile(exp, f)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    parser.add_argument('--real', action="store_true", default=False)
    parser.add_argument('--host', action="store_true",
                        default="https://test.encodedcc.org")
    parser.add_argument('--keyFnp', action="store_true",
                        default=os.path.expanduser('~/.encode.txt'))
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.real:
        args.host = "https://www.encodeproject.org"

    s = SubmitHiCLiftOver(args)

    if 1:
        return s.runLiftover()
    if 0:
        return s.runSubmit()


if __name__ == "__main__":
    sys.exit(main())
