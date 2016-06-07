#!/usr/bin/env python

import os, sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../metadata/utils'))
from utils import Utils, cat, numLines, printWroteNumLines
from files_and_paths import Dirs, Tools, AllHumanDataset

inFnp = sys.argv[1]
print inFnp

genePredFnp = os.path.splitext(inFnp)[0] + '.genePred'
cmds = [Dirs.ToolsFnp("ucsc.2016-02Feb-16/gtfToGenePred"),
        inFnp, genePredFnp]
Utils.runCmds(cmds, True)
print genePredFnp

bedPlusFnp = os.path.splitext(inFnp)[0] + '.bedPlus'
cmds = [Dirs.ToolsFnp("ucsc.2016-02Feb-16/genePredToBigGenePred"),
        genePredFnp, bedPlusFnp]
Utils.runCmds(cmds, True)
print bedPlusFnp

bedPlusSortedFnp = os.path.splitext(inFnp)[0] + '.sorted.bedPlus'
cmds = [cat(bedPlusFnp),
        '|', "LC_COLLATE=C sort -k1,1 -k2,2n",
        '>', bedPlusSortedFnp]
Utils.runCmds(cmds, True)
print bedPlusFnp

bigPredFnp = os.path.splitext(inFnp)[0] + '.bigGenePred'
cmds = [Dirs.ToolsFnp("ucsc.2016-02Feb-16/bedToBigBed"),
        "-type=bed12+8",
        "-tab",
        "-as=" + Dirs.ToolsFnp("ucsc.2016-02Feb-16/bigGenePred.as"),
        bedPlusSortedFnp,
        AllHumanDataset.chr_lengths,
        bigPredFnp]
Utils.runCmds(cmds, True)
print bigPredFnp

print("done", inFnp)
