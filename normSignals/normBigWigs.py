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

def process(idx, job):
    args = job[0]
    exp = job[1]
    try:
        bigWigFnp, bwAssembly = exp.getSingleBigWigSingleFnp(args)
        if not bigWigFnp:
            print exp.getSingleBigWigSingleFnp(args)
            print "missing", exp
        else:
            cmds = [normBin,
                    "--assembly=" + bwAssembly,
                    bigWigFnp]
            return Utils.runCmds(cmds)
    except Exception, e:
        return "bad " + str(e)

def build(args):
    jr = JobRunner(scriptFnp = os.path.realpath(__file__),
                   jobType = PythonJob,
                   cpus = args.j)

    epigenomes = WebEpigenomesLoader(args)
    for assembly in ["hg19", "mm10", "mm9"]:
        for assays in ["H3K27ac", "DNase"]:
            epis = epigenomes.GetByAssemblyAndAssays(assembly, assays)
            for epi in epis.epis:
                for exp in epi.exps():
                    jr.append([[args, exp]])

    if args.test:
        return jr.runOne(process)

    if args.local:
        return jr.run(runJob)

    jobOptions = {"mem" : 64000,
                  "time" : "3:59",
                  "cores" : 1,
                  "queue" : "short" }

    jr.cluster("/project/umw_zhiping_weng/encyc/norm", jobOptions)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--process', action="store_true", default=True)
    parser.add_argument('--local', action="store_true", default=False)
    parser.add_argument('--rebuild', action="store_true", default=False)
    parser.add_argument('--test', action="store_true", default=False)
    parser.add_argument('--job', type=int, default=0)
    parser.add_argument('-j', type=int, default=4)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    return build(args)

if __name__ == '__main__':
    main()
