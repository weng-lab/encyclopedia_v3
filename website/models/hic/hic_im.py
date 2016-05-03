#!/usr/bin/env python

import sys, os
import h5py
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../metadata/utils'))
from utils import Utils

class HiCInteractionMatrix:
    def __init__(self):
        self.params = []

    def makeImg(fnp, chrom, start, numBins, res):
        fnp = "/home/mjp/hic-Dekker/hdf5/a549/ENCODE3-A549-HindIII__hg19__hdf/C-40000/iced/ENCODE3-A549-HindIII__hg19__genome__C-40000-iced.hdf5"

        f = h5py.File(fnp, "r")

        #print f.keys()
        chrs = dict((c, i) for i, c in enumerate(f["chrs"]))
        chromIdx = chrs[chrom]
        chr_bin_range = f["chr_bin_range"][chromIdx]
        print chr_bin_range

        binIdx = start / res
        idx = chr_bin_range[0] + binIdx
        print chr_bin_range, idx, numBins

        im = f["interactions"]
        m = im[idx : idx + numBins, idx : idx + numBins]
        print m.shape

        pngFnp = fnp + ".png"
        saveMatrix(pngFnp, np.matrix(m))
        rotate(pngFnp)

    def saveMatrix(fnp, im):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        m = im.clip(0, 1)
        m = np.triu(m)
        m[m == 0.0] = np.nan

        plt.figure(1)
        plt.imshow(m, interpolation='nearest')
        plt.axis('off')
        plt.savefig(fnp)

    def rotate(fnp):
        from PIL import Image

        s = Image.open(fnp)
        s = s.convert('RGBA')
        s = s.rotate(45, expand=1)

        # http://stackoverflow.com/a/5253554
        # a white image same size as rotated image
        fff = Image.new('RGBA', s.size, (255,)*4)

        # create a composite image using the alpha layer of rot as a mask
        out = Image.composite(s, fff, s)

        # save your work (converting back to mode='1' or whatever..)
        out.convert(s.mode).save(fnp)

        # trim
        cmds = ["convert",
                fnp,
                "-trim +repage",
                fnp + ".trim.png"]
        Utils.runCmds(cmds)

def main():
    im = HiCInteractionMatrix()
    im.makeImg("", "chr1", 40000, 2002, 40000)

if __name__ == '__main__':
    main()


