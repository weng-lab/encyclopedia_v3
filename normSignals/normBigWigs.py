#!/usr/bin/env python

import os, sys, json, argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '../website/'))
from web_epigenomes import WebEpigenomesLoader
from helpers_trackhub import bigWigFilters

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils, printWroteNumLines
from job_runner import JobRunner, PythonJob
from metadataws import MetadataWS
from files_and_paths import Datasets

normBin = os.path.join(os.path.realpath(os.path.dirname(__file__)), "bin/normBigWig")
if not os.path.exists(normBin):
    print "missing", normBin
    sys.exit(1)

def process(args, expID):
    exp = MetadataWS.exp(expID)
    try:
        bigWigs = bigWigFilters(args.assembly, exp.files)
        if not bigWigs:
            return
        bigWig = bigWigs[0]
        bigWigFnp = bigWig.fnp()
        if os.path.exists(bigWig.normFnp()):
            print "skipping", exp
            return 0
        else:
            print "missing", bigWig.normFnp()
        bwAssembly = bigWig.assembly
        if not bigWigFnp:
            print exp.getSingleBigWigSingleFnp(args)
            print "missing", exp
        else:
            cmds = [normBin,
                    "--assembly=" + bwAssembly,
                    "--bwFnp=" + bigWig.normFnp(),
                    bigWigFnp]
            print " ".join(cmds)
            print Utils.runCmds(cmds)
            return 0
    except Exception, e:
        print "bad " + str(e)
    return 1

def build(args):
    jr = JobRunner(cpus = args.j)

    epigenomes = WebEpigenomesLoader(args)
    for assembly in ["hg19", "mm10"]:
        for assays in ["H3K27ac", "DNase"]:
            epis = epigenomes.GetByAssemblyAndAssays(assembly, assays)
            for epi in epis.epis:
                for exp in epi.exps():
                    if exp.encodeID.startswith("EN"):
                        #print exp.encodeID
                        jr.append([
                                [os.path.realpath(__file__), "--job", exp.encodeID,
                                 "--assembly", assembly, "--process"]])
                    else:
                        # ROADMAP
                        bigWig = exp.files[0]
                        if not bigWig:
                            print "missing", exp
                        else:
                            if not os.path.exists(bigWig.normFnp()):
                                jr.append([
                                        [normBin,
                                         "--assembly=" + bigWig.assembly,
                                         "--bwFnp=" + bigWig.normFnp(),
                                         bigWig.fnp()]])
    if args.test:
        return jr.runOne()

    if args.local:
        return jr.run()

    jobOptions = {"mem" : 64000,
                  "time" : "3:59",
                  "cores" : 2,
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
    parser.add_argument('--assembly', type=str, default="hg19")
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
