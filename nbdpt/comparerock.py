"""
Compare a mass function from rockstar to theoretical mass functions
from Behroozi+ 

import nbdpt.comparerock as rc
rc.compare()

keyword arguments
------
  nrows    - if doing multiple subplots, the number of rows
  ncolumns - if doing multiple subplots, the number of columns
  nbins    - number of bins for histogram of halo masses
  filename - filename to save plots to
  subplot  - None
  z        - a list of the redshifts of comparison

"""
from . import amiga 
from . import rockstar
from . import nptipsyreader
import numpy as np
import glob
import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb

def readoutputs(tipsyfilename):
    # read in tipsy file and corresponding rockstar output
    output = tipsyfilename.split('.')[-1]
    f = open('snaps.txt')
    snaps = [snaps.strip('\n') for snaps in f]
    f.close()
    snaps = np.array(snaps)
    try: 
        outnumber = np.where(snaps == output)[0][0]
        rock = rockstar.readout.read('out_' + str(outnumber) + '.list')
        return rock
    except IndexError:
        print "Tipsyfile not Rockstar'd"
        return None

def readmf(tipsyfilename, h):
    # read in the theoretical mass function
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    if np.isclose(h,0.6777,rtol=1e-3): mfname = 'data/mf_planck13.dat'
    if np.isclose(h,0.73,rtol=1e-3):   mfname = 'data/mf_wmap3.dat'
    if np.isclose(h,0.70,rtol=1e-3):   mfname = 'data/mf_wmap7.dat'
    mfile = os.path.join(_ROOT, mfname)
    f = open(mfile)
    names = f.readline().strip('\n#').split()
    f.close()
    mf = np.genfromtxt(mfile, names=names, skiprows=1)
    return mf


def plot(rock, mf, h, boxsize, thiscolor, redshift, roundedz, thisplot, nbins):
    #plot the rockstar output and theoretical mass function 
    rocknz = rock['M200c'] > 0.
    logmass = np.log10(rock[rocknz]['M200c']/h)
    hist, bins = np.histogram(logmass, bins=nbins)
    normalization = (bins[1:] - bins[:-1])*(boxsize)**3.
    err = 1. + (hist + 0.75)**0.5
    thisplot.errorbar(10.**((bins[1:]+bins[:-1])/2.), 
                      hist/normalization, yerr=err/normalization, 
                      linewidth=2, label='z~' + '{:.2f}'.format(redshift)+' Rockstar', 
                      color=thiscolor, drawstyle='steps-mid')
    mfm  = mf['M']/h
    dndm = mf['dndlog10m_z' + str(int(roundedz)) + '0']*h**3.
    thisplot.plot(mfm, dndm, label='Behroozi+ z~' + str(int(roundedz)), linewidth=2, color=thiscolor)

def compare(nrow=1, ncolumns=1, filename='massfunctions.png', subplot=None, 
            z = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], nbins=10.):    
    #compare rockstar and theoretical
    fontsize = 10
    threshold = 0.1 #allowed deviation from theoretical mf redshift for comparison
    tipsyfiles = [f for f in os.listdir('.') if re.match('^[\d]*$', f.split('.')[-1])]
    tipsyfiles.sort()
    tipsyfiles = tipsyfiles[::-1]
    colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(tipsyfiles))))
    for tipsyfilename in tipsyfiles:
        sim = nptipsyreader.Tipsy(tipsyfilename)
        sim._read_param()
        (t, nbodies, ndim, ngas, dark, nstar) = sim.unpack_header()
        redshift = 1./t-1
        boxsize = np.float(sim.paramfile['dKpcUnit'])/1e3
        roundedz = round(redshift/threshold)*threshold
        if roundedz in z:
            if subplot: thisplot=subplot
            else: fig, thisplot = plt.subplots()
            rock = readoutputs(tipsyfilename)
            if rock is None: continue 
            mf = readmf(tipsyfilename, sim.h)
            thiscolor = next(colors)
            plot(rock, mf, sim.h, boxsize, thiscolor, redshift, roundedz, thisplot, nbins)    
            thisplot.set_xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
            thisplot.set_ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$', fontsize=fontsize)
            thisplot.set_xscale('log')
            thisplot.set_yscale('log')
            thisplot.legend()
            plt.show()
            #pdb.set_trace()
    if not(subplot): plt.savefig(filename)
    plt.show()
