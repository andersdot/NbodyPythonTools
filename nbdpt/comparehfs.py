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
    output = tipsyfilename.split('.')[-1]
    
    f = open('snaps.txt')
    snaps = [snaps.strip('\n') for snaps in f]
    f.close()
    snaps = np.array(snaps)
    
    ahf = amiga.newAhfGrpStat.AHF(tipsyfilename)
    ahfhalos = ahf.halos()
    outnumber = np.where(snaps == output)[0][0]
    rock = rockstar.readout.read('out_' + str(outnumber) + '.list')
    return ahfhalos, rock

def readmf(tipsyfilename, h):
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


def plot(ahf, rock, mf, h, boxsize, redshift, roundedz, thisplot):
    fig, ax = plt.subplots()
    rocknz = rock['M200c'] > 0.
    amiganz = ahf['Mvir4'] > 0.
    logrockmass = np.log10(rock['M200c'][rocknz]/h)
    logamigamass = np.log10(ahf['Mvir4'][amiganz]/h)
    
    hist, bins = np.histogram(logamigamass, bins=5)
    normalization = (bins[1:] - bins[:-1])*(boxsize)**3.
    err = 1. + (hist + 0.75)**0.5
    ax.errorbar(10.**((bins[1:]+bins[:-1])/2.),
                hist/normalization, yerr=err/normalization,
                linewidth=2, label='Amiga', color='grey')
    
    am = thisplot.errorbar(10.**((bins[1:]+bins[:-1])/2.), 
                           hist/normalization, yerr=err/normalization, 
                           linewidth=2, label='Amiga', color='grey')
    
    hist, bins = np.histogram(logrockmass, bins=5)
    normalization = (bins[1:] - bins[:-1])*(boxsize)**3.
    err = 1. + (hist + 0.75)**0.5
    ro = thisplot.errorbar(10.**((bins[1:]+bins[:-1])/2.),
                           hist/normalization, yerr=err/normalization,
                           linewidth=2, label='Rockstar', color='blue')
    ax.errorbar(10.**((bins[1:]+bins[:-1])/2.),
                hist/normalization, yerr=err/normalization,
                linewidth=2, label='Rockstar', color='blue')
    
    mfm  = mf['M']/h
    fixedz = int(roundedz*10.)
    if fixedz < 10: dndm = mf['dndlog10m_z0' + str(int(roundedz*10.))]*h**3.
    else: dndm = mf['dndlog10m_z' + str(int(roundedz*10.))]*h**3.
    ax.plot(mfm, dndm, label='Berhoozi+', linewidth=2, color='k')
    mf = thisplot.plot(mfm, dndm, label='Berhoozi+', linewidth=2, color='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fontsize=12
    ax.set_xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
    ax.set_ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$',
                        fontsize=fontsize)
    ax.set_title('z~' + str(int(roundedz)) + ' z='+'{:.2f}'.format(redshift),fontsize=10)
    ax.legend(prop={'size':10})
    plt.show()
    return am, ro,mf
def compareall():    

    fontsize = 10
    threshold = 0.1 #allowed deviation from theoretical mf redshift for comparison
    tipsyfiles = [f for f in os.listdir('.') if re.match('^[\d]*$', f.split('.')[-1])]
    tipsyfiles.sort()
    colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(tipsyfiles))))
    sub = 0
    jj = 0
    fig, ax = plt.subplots(4,4)
    fig.tight_layout()

    for tipsyfilename in tipsyfiles:
        sim = nptipsyreader.Tipsy(tipsyfilename)
        sim._read_param()
        (t, nbodies, ndim, ngas, dark, nstar) = sim.unpack_header()
        redshift = 1./t-1
        boxsize = np.float(sim.paramfile['dKpcUnit'])/1e3
        roundedz = round(redshift/threshold)*threshold
        print roundedz, jj
        z = [1.0] #[0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        if (roundedz in z) and (jj < 16):
            ahf, rock = readoutputs(tipsyfilename)
            mf = readmf(tipsyfilename, sim.h)
            thisplot = ax[sub/4, sub%4]
            print sub, sub/4, sub%4
            #thisplot.xaxis.set_major_formatter(plt.LogLocator(subs=2.0, numticks=4))
            #thisplot.yaxis.set_major_formatter(plt.LogLocator(subs=2.0, numticks=4))
            am, ro, mf = plot(ahf, rock, mf, sim.h, boxsize, redshift, roundedz, thisplot)    
            
            thisplot.set_xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
            thisplot.set_ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$', 
                                fontsize=fontsize)
            thisplot.set_xscale('log')
            thisplot.set_yscale('log')
            thisplot.set_title('z~' + str(int(roundedz)) + ' z='+'{:.2f}'.format(redshift),fontsize=10)
            sub += 1
            jj += 1
            pdb.set_trace()
            plt.show()
    fig.legend((am, ro, mf), ('Amiga', 'Rockstar', 'Behroozi+'), 1)
    plt.savefig('cosmo50p.512g.massfunctions.png')
    plt.show()
