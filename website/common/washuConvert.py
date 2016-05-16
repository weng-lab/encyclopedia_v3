#!/usr/bin/env python

import os, sys, shutil
from joblib import Parallel, delayed

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils'))
from utils import Utils, printWroteNumLines
from files_and_paths import Dirs

def job(args):
    fn = args[0]
    inFnp = args[1]
    outD = args[2]

    outFnp = os.path.join(outD, fn)
    with open(inFnp) as inF:
        with open(outFnp, 'w') as outF:
            for idx, line in enumerate(inF):
                toks = line.rstrip().split('\t')
                attrs = "id:"+str(idx)+',name:"'+toks[3]+'"'
                if 8 == len(toks):
                    attrs += ",struct:{{thick:[[{s},{e}],],}}".format(s=toks[6], e=toks[7])
                out = toks[:3] + [attrs]
                outF.write("\t".join(out) + '\n')
    Utils.sortFile(outFnp)
    printWroteNumLines(outFnp)

    cmds = ["bgzip", '-f', outFnp]
    Utils.runCmds(cmds)

    cmds = ["tabix", '-f', outFnp + '.gz']
    Utils.runCmds(cmds)

    print("wrote", inFnp, outFnp)

def runDir(d):
    outD = os.path.join(d, "washu")
    Utils.mkdir_p(outD)

    jobs = []

    for fn in os.listdir(d):
        if not fn.endswith(".bed"):
            continue
        inFnp = os.path.join(d, fn)
        jobs.append([fn, inFnp, outD])

    ret = Parallel(n_jobs = 4)(delayed(job)(j)
                               for idx, j in enumerate(jobs))

def main():
    runDir(Dirs.promoterTracks)
    runDir(Dirs.enhancerTracks)


if __name__ == '__main__':
    main()
