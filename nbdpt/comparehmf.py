"""    
import nbdpt.comparehmf as chmf
chmf.compareall('cosmo12PLK.256.000972', 'cosmo8.33PLK.256g2bwK1C52.000970', amiga=True)

output: massfunctions.png

"""

from . import nchiladareader as nch
from . import readstat as rds
from . import nptipsyreader as npt
from . import findtipsy as ft
import numpy as np
import glob
import os
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb
from astroML.plotting import hist as histogram 

fontsize=12

def readstat(simfilename, amiga=False):
    
    if amiga: stat = rds.readstat(simfilename + '.amiga.stat', amiga=amiga)
    else: stat = rds.readstat(simfilename + '.rockstar.stat', amiga=amiga)
    return stat

def readmf(simfilename):
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    nchil = nch.Nchilada(simfilename)
    nchil.read_param()
    h = nchil.h
    print 'h = ' + str(h)
    print 'Filename: ' + simfilename + '   Paramfilename: ' + nchil.paramfilename  
    if np.isclose(h,0.6777,rtol=1e-3): mfname = 'data/mf_planck13.dat'
    if np.isclose(h,0.73,rtol=1e-3):   mfname = 'data/mf_wmap3.dat'
    if np.isclose(h,0.70,rtol=1e-3):   mfname = 'data/mf_wmap7.dat'
    print 'Assuming ' + mfname.split('.')[0].split('/')[-1].split('_')[-1] + ' cosmology'
    mfile = os.path.join(_ROOT, mfname)
    f = open(mfile)
    names = f.readline().strip('\n#').split()
    f.close()
    mf = np.genfromtxt(mfile, names=names, skiprows=1)
    return mf


def plotstat(simfilename, thisplot, amiga=False, bins=False):
    
    stat = readstat(simfilename, amiga=amiga)

    r = re.compile('cosmo(.*?)PLK')
    m = r.search(simfilename)
    if m: boxsize = float(m.group(1))
    else: 
        r = re.compile('cosmo(.*?)cmb')
        m = r.search(filename)
        boxsize = float(m.group(1))
    if boxsize: print 'Boxsize is ' + str(boxsize) + ' Mpc'
    else: print "Couldn't extract boxsize from simulation filename"

    try:
        nchil = nch.Nchilada(simfilename)
        nchil.read_param()
        h = nchil.h
        a, nbodies, ndim, code = nchil.unpack_header('dark', 'mass')
        redshift = 1./a -1.
    except IOError:
        try:
            tipsy = npt.Tipsy(simfilename)
            tipsy._read_param()
            h = tipsy.h
            a, nbodies, ndim, nsph, ndark, nstar = tipsy.unpack_header()
            redshift = 1./a -1.
        except IOError:
            redshift = input("Can't find " + simfilename + ", what is the redshift?")

    print 'h = ' + str(h)
    print 'The redshift is ' + str(redshift)
    statnz = (stat['mvir'] > 0.) & (stat['npart'] > 512.)

    if amiga: 
        logmass = np.log10(stat['mvir'][statnz])
        label = 'Amiga'
    else: 
        logmass = np.log10(stat['mvir'][statnz]/h)
        label = 'Rockstar'

    #hist, bins = np.histogram(logmass, bins=5)
    dumbf, dumb = plt.subplots(1)
    if bins is not None:
        hist, bins = np.histogram(logmass, bins=bins)
    else:
        output = histogram(logmass, bins='blocks', ax=dumb)
        hist = output[0]
        bins = output[1]
    normalization = (bins[1:] - bins[:-1])*(boxsize)**3.
    err = 1. + (hist + 0.75)**0.5
    am = thisplot.errorbar(10.**((bins[1:]+bins[:-1])/2.), 
                           hist/normalization, yerr=err/normalization, 
                           linewidth=2, label=simfilename.split('.')[0])
    thisplot.set_xscale('log')
    thisplot.set_yscale('log')
    thisplot.set_xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
    thisplot.set_ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$',
                        fontsize=fontsize)
    thisplot.set_title('z='+'{:.2f}'.format(redshift),fontsize=fontsize)
    thisplot.legend(prop={'size':10})

    return thisplot, bins

def plotmf(simfilename, thisplot):

    try:
        nchil = nch.Nchilada(simfilename)
        nchil.read_param()
        h = nchil.h
        a, nbodies, ndim, code = nchil.unpack_header('dark', 'mass')
        redshift = 1./a -1.
    except IOError:
        try:
            tipsy = npt.Tipsy(simfilename)
            tipsy._read_param()
            h = tipsy.h
            a, nbodies, ndim, nsph, ndark, nstar = tipsy.unpack_header()
            redshift = 1./a -1.
        except IOError:
            redshift = input("Can't find " + simfilename+ ", what is the redshift?")
    print 'h = ' + str(h)
    print 'The redshift is ' + str(redshift)
    threshold = 0.1
    roundedz = round(redshift/threshold)*threshold
    
    mf = readmf(simfilename)
    mfm  = mf['M']/h
    fixedz = int(roundedz*10.)
    if fixedz < 10: dndm = mf['dndlog10m_z0' + str(int(roundedz*10.))]*h**3.
    else: dndm = mf['dndlog10m_z' + str(int(roundedz*10.))]*h**3.
    mf = thisplot.plot(mfm, dndm, label='Berhoozi+', linewidth=2, color='k')
    thisplot.set_xscale('log')
    thisplot.set_yscale('log')
    thisplot.set_xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
    thisplot.set_ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$',
                        fontsize=fontsize)
    thisplot.set_title('z~' + str(int(roundedz)),fontsize=fontsize)
    thisplot.legend(prop={'size':10})

    return thisplot

def compareall(simfilenames=False, amiga=False):    

    if not(simfilenames): simfilenames = ft.find()
    colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(simfilenames))))
    fig, ax = plt.subplots(1)
    #fig.tight_layout()
    bins = None
    for sim in simfilenames:
        print sim, amiga
        #z = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

        #if (roundedz in z):
        ax, bins = plotstat(sim, ax, amiga=amiga, bins=bins)
        ax = plotmf(sim, ax)
    fig.savefig('massfunctions.png')


