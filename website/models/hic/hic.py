#!/usr/bin/env python

import h5py
import numpy as np

fnp = "/home/mjp/hic-Dekker/hdf5/a549/ENCODE3-A549-HindIII__hg19__hdf/C-10000000/iced/ENCODE3-A549-HindIII__hg19__genome__C-10000000-iced.hdf5"

f = h5py.File(fnp, "r")

print f.keys()
chrs = f["chrs"]
chr_bin_range = f["chr_bin_range"]

im = f["interactions"]

print im.shape

for idx, c in enumerate(chrs):
    cbr = chr_bin_range[idx]
    #print c, ",".join([str(x) for x in im[cbr[0] : cbr[1]]])

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
m = np.matrix(im).clip(0, 1)

plt.figure(1)
plt.imshow(m, interpolation='nearest')
plt.savefig(fnp + ".png")
