#!/usr/bin/env python

import os, sys, json, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../website/'))
from web_epigenomes import WebEpigenomesLoader

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils, printWroteNumLines
from job_runner import JobRunner, PythonJob
from metadataws import MetadataWS
from files_and_paths import Datasets

normBin = "/home/purcarom/annotations/normSignals/bin/normBigWig"

def process(args, expID):
    exp = MetadataWS.exp(expID)
    print exp
    try:
        bigWigFnp, bwAssembly = exp.getSingleBigWigSingleFnp(args)
        if not bigWigFnp:
            print exp.getSingleBigWigSingleFnp(args)
            print "missing", exp
        else:
            cmds = [normBin,
                    "--assembly=" + bwAssembly,
                    bigWigFnp]
            print "running", cmds
            print Utils.runCmds(cmds)
            return 0
    except Exception, e:
        print "bad " + str(e)
        return 1
    return 0

def build(args):
    jr = JobRunner(cpus = 1)

    epigenomes = WebEpigenomesLoader(args)
    for assembly in ["hg19", "mm10", "mm9"]:
        for assays in ["H3K27ac", "DNase"]:
            epis = epigenomes.GetByAssemblyAndAssays(assembly, assays)
            for epi in epis.epis:
                for exp in epi.exps():
                    if exp.encodeID.startwith("EN"):
                        jr.append([
                                [__file__, "--job", exp.encodeID,
                                 "--process"]])
                    else:
                        # ROADMAP
                        bigWigFnp, bwAssembly = exp.getSingleBigWigSingleFnp()
                        if not bigWigFnp:
                            print "missing", exp
                        else:
                            cmds = [normBin,
                                    "--assembly=" + bwAssembly,
                                    bigWigFnp]

    if args.test:
        return jr.runOne()

    if args.local:
        return jr.run()

    jobOptions = {"mem" : 64000,
                  "time" : "3:59",
                  "cores" : 1,
                  "queue" : "short" }

    jr.cluster("/project/umw_zhiping_weng/encyc/norm", jobOptions)

def runJob(args):
    return process(args, args.job)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--process', action="store_true", default=True)
    parser.add_argument('--local', action="store_true", default=False)
    parser.add_argument('--rebuild', action="store_true", default=False)
    parser.add_argument('--test', action="store_true", default=False)
    parser.add_argument('--job', type=str, default="")
    parser.add_argument('-j', type=int, default=4)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    if args.job:
        return runJob(args)

    return build(args)

if __name__ == '__main__':
    sys.exit(main())
