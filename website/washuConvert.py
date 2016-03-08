#!/usr/bin/env python

import os, sys, shutil

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils'))
from utils import Utils
from files_and_paths import Dirs

def main():
    d = os.path.join("/home/purcarom/public_html/annotations_demo",
                     "Enhancer-Prediction-Tracks")
    outD = os.path.join(d, "washu")
    Utils.mkdir_p(outD)

    for fn in os.listdir(d):
        if not fn.endswith(".bed"):
            continue
        inFnp = os.path.join(d, fn)
        print(inFnp)

        outFnp = os.path.join(outD, fn)
        shutil.copyfile(inFnp, outFnp)

        cmds = ["bgzip", '-f', outFnp]
        Utils.runCmds(cmds)

        cmds = ["tabix", '-f', outFnp + '.gz']
        Utils.runCmds(cmds)

if __name__ == '__main__':
    main()
