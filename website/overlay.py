#!/usr/bin/env python

import os, sys, json, argparse
from collections import defaultdict
from itertools import groupby
import cPickle as pickle

import numpy as np

from natsort import natsorted

from roadmap import RoadmapMetadata

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from metadataws import MetadataWS
from files_and_paths import Datasets
from web_epigenomes import WebEpigenomesLoader

from helpers_trackhub import Track, PredictionTrack, BigGenePredTrack, BigWigTrack, officialVistaTrack, bigWigFilters
from trackhub import TempWrap

def process(assembly, exp):
    bigWigs = bigWigFilters(assembly, exp.files)
    if not bigWigs:
        if "mm10" == assembly:
            bigWigs = [TempWrap(exp.encodeID,
                                eidToBigWigFileID[exp.encodeID])]
    if not bigWigs:
        print "missing bigWigs for " + exp.encodeID
        return

    if 0 and len(bigWigs) > 1:
        print "weird", bigWigs
        sys.exit(1)

    bigWig = bigWigs[0]
    if "eid" in bigWig.fnp():
        return # ignore Roadmap for now
    print bigWig.fnp()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action="store_true", default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    we = WebEpigenomesLoader(args)

    for assembly in ["mm9", "mm10", "hg19"]:
        for assays in ["H3K27ac", "DNase"]:
            wepis = we.GetByAssemblyAndAssays(assembly, assays)
            for epi in wepis.epis:
                for exp in epi.exps():
                    process(assembly, exp)

if __name__ == '__main__':
    main()

