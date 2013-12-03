import os
import numpy as np
import pylab
from .. import nptipsyreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
from .. import statfile
import pdb
#######################
## Compare Redshifts ##
#######################

def plot():
    fontsize = 18
    linewidth = 2
    
    ylim = (-6, 2.5)
    xlim = (6, 15.5)
    z = [0, 0.1, 0.5, 1, 2, 3, 4] #, 5, 6, 7, 8, 9, 10]
    
    prefix = '.'.join(glob.glob('*.iord')[0].split('.')[:-2])
    snapfile = 'snaps.txt'
    f = open(snapfile)
    snaps = [snap.strip('\n') for snap in f]
    f.close()
    testoutput = prefix + '.' + str(snaps[0])
    testsim = nptipsyreader.Tipsy(testoutput)
    testsim._read_param()

    _ROOT = os.path.abspath(os.path.dirname(__file__))
    
    if np.isclose(testsim.h,0.6777,rtol=1e-3): mfname = 'data/mf_planck13.dat'
    if np.isclose(testsim.h,0.73,rtol=1e-3):   mfname = 'data/mf_wmap3.dat'
    if np.isclose(testsim.h,0.70,rtol=1e-3):   mfname = 'data/mf_wmap7.dat'
    
    mfile = os.path.join(_ROOT, mfname)

    f = open(mfile)
    names = f.readline().strip('\n#').split()
    f.close()
    mf = np.genfromtxt(mfile, names=names, skiprows=1)
    
    rockoutputs = glob.glob('out_*.list')
    rockoutputs.sort()
    colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(rockoutputs))))
    
    for i in range(len(rockoutputs)):
        output = prefix +'.'+ str(snaps[i]) 
        tipsy = nptipsyreader.Tipsy(output)
        tipsy._read_param()
        h = tipsy.h
        box = np.float(tipsy.paramfile['dKpcUnit'])/1e3
        (t, nbodies, ndim, nsph, ndark, nstar) = tipsy.unpack_header()
        redshift = 1./t -1.
        roundedz = np.int(np.round(redshift))
        if roundedz in z:
            f = open(rockoutputs[i])
            names = f.readline().strip('\n#').split(' ')
            dtype = np.concatenate((['int', 'int'], (36*'float ').split(' ')[:-1]))
            f.close()
            stat = np.genfromtxt(rockoutputs[i], dtype=dtype, names=names, skiprows=1)
            nonzero = stat['M200c'] > 0.
            logmass = np.log10(stat['M200c'][nonzero]/h)
            hist, bins = np.histogram(logmass, bins=5)
            normalization = (bins[1:] - bins[:-1])*(box)**3.
            err = 1. + (hist + 0.75)**0.5
            thiscolor = next(colors)
            plt.errorbar(10.**((bins[1:]+bins[:-1])/2.), 
                         hist/normalization, yerr=err/normalization, 
                         linewidth=2, label='z~' + '{:.2f}'.format(redshift), color=thiscolor)
            
            
            mfm  = mf['M']/h
            dndm = mf['dndlog10m_z' + str(roundedz) + '0']*h**3.
            plt.plot(mfm, dndm, label='Behroozi+ z~' + str(roundedz), linewidth=2, color=thiscolor)
            
    plt.xlabel('$\mathrm{log_{10} M_{tot}}$', fontsize=fontsize)
    plt.ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$', fontsize=fontsize)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig('cosmo50p.512.massfunction.png')

if __name__=='__main__':plot()
