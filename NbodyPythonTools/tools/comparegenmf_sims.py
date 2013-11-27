import numpy as np
import pylab
import nptipsyreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import statfile
#######################
## Compare Redshifts ##
#######################

fontsize = 18
linewidth = 2

ylim = (-6, 2.5)
xlim = (6, 15.5)
z = [0, 1, 2, 3, 4, 5]
mfdir = '/u/sciteam/lmanders/data/'
mfname = 'plk13.allz.mf'
f = open(mfdir+mfname)
names = f.readline().strip('\n#').split()
f.close()
mf = np.genfromtxt(mfdir + mfname, names=names, skiprows=1)

outputs = glob.glob('*.stat')
outputs.sort()
colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(outputs))))

for i in range(len(outputs)):
    file = ('.').join(outputs[i].split('.')[:-4])
    tipsy = nptipsyreader.Tipsy(file)
    tipsy._read_param()
    h = tipsy.h
    box = np.float(tipsy.paramfile['dKpcUnit'])/1e3
    (t, nbodies, ndim, nsph, ndark, nstar) = tipsy.unpack_header()
    redshift = 1./t -1.
    roundedz = np.int(np.round(redshift))
    if roundedz in z:
        statfilename = glob.glob(file + '.fof.*.stat')[0]
        print statfilename
        stat = statfile.readfof(statfilename)
        nonzero = stat['mtot'] > 0.
        logmass = np.log10(stat['mtot'][nonzero])
        hist, bins = np.histogram(logmass, bins=5)
        normalization = (bins[1:] - bins[:-1])*(box)**3.
        err = 1. + (hist + 0.75)**0.5
        thiscolor = next(colors)
        plt.errorbar(10.**((bins[1:]+bins[:-1])/2.), 
                     hist/normalization, yerr=err/normalization, 
                     linewidth=2, label='z~' + str(roundedz), color=thiscolor)

        
        mfm  = mf['M']/h
        dndm = mf['dndlog10m_z' + str(roundedz) + '0']*h**3.
        plt.plot(mfm, dndm, label='Behroozi+ z~' + str(roundedz), linewidth=2, color=thiscolor)
        
plt.xlabel('$\mathrm{log_{10} M_{tot}}$)', fontsize=fontsize)
plt.ylabel('$\Phi \; [\mathrm{Mpc}^{-3} \; \mathrm{log_{10} M_{tot}}^{-1}]$', fontsize=fontsize)
plt.xscale('log')
plt.yscale('log')
plt.show()
