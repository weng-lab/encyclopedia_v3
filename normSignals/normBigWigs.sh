#!/usr/bin/env python

import os, sys, json, argparse

from web_epigenomes import WebEpigenomesLoader

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils, printWroteNumLines
from metadataws import MetadataWS
from files_and_paths import Datasets

def build(args):
    fnps = []

    epigenomes = WebEpigenomesLoader(args)
    for assembly in ["hg19", "mm10", "mm9"]:
        for assays in ["H3K27ac", "DNase"]:
            epis = epigenomes.GetByAssemblyAndAssays(assembly, assays)
            for epi in epis.epis:
                for exp in epi.exps():
                    try:
                        bigWigFnp, bedAssembly = exp.getSingleBigWigSingleFnp(args)
                        if not bigWigFnp:
                            print exp.getSingleBigWigSingleFnp(args)
                            print "missing", exp
                        else:
                            fnps.append(bigWigFnp)
                    except Exception, e:
                        print str(e)
                        print "bad exp:", exp

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--process', action="store_true", default=True)
    parser.add_argument('--local', action="store_true", default=False)
    parser.add_argument('--rebuild', action="store_true", default=False)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    return build(args)

if __name__ == '__main__':
    main()
