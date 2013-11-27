import readstat
import nptipsyreader
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb

linewidth = 3.
gonzredshift = [4.0, 5.0, 6.0, 7.0] 
gonzdndm = [np.array((-1.90, -1.96, -2.30, -2.53, -3.14, -3.80, -4.43)),
             np.array((-2.21, -2.27, -2.60, -2.84, -3.46, -4.12, -4.81)),
             np.array((-2.09, -2.16, -2.55, -2.81, -3.52, -4.22, -4.97)),
             np.array((-2.15, -2.23, -2.70, -3.01, -3.80, -4.53, -1e6))]
 
gonzmasses = np.array((7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75))
gonzbinwidth = np.array((0.5, 0.5,  0.5,  0.5,  0.5,   0.5,   0.5))

czakredshift = [0.5, 0.75, 1., 1.25, 1.5, 2.0]
czakdndm = [np.array((-1.53, -1.60, -1.76, -1.86, -2.00, -2.12, -2.21, -2.25, -2.35, -2.45, -2.55, -2.82, -3.32, -1e6)),
            np.array((-1e6 , -1.70, -1.86, -2.01, -2.10, -2.23, -2.39, -2.45, -2.45, -2.52, -2.59, -2.93, -3.47, -1e6 )),
            np.array((-1e6 , -1e6 , -1.99, -2.14, -2.24, -2.29, -2.48, -2.59, -2.73, -2.64, -2.72, -3.01, -3.62, -1e6 )),
            np.array((-1e6 , -1e6 , -2.02, -2.14, -2.28, -2.46, -2.53, -2.61, -2.68, -2.71, -2.84, -3.12, -3.65, -4.99)),
            np.array((-1e6 , -1e6 , -1e6 , -2.20, -2.31, -2.41, -2.54, -2.67, -2.76, -2.87, -3.03, -3.13, -3.56, -4.27)),
            np.array((-1e6 , -1e6 , -1e6 , -1e6 , -2.53, -2.50, -2.63, -2.74, -2.91, -3.07, -3.35, -3.54, -3.89, -4.41))]
czakmasses =  np.array((8.25, 8.50, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.50))

nbins = [np.array(( 5.423,  6.564,  7.704)),
        np.array(( 5.404,  6.473,  7.542,  8.611)),
        np.array(( 4.582,  5.833,  7.085,  8.336,   9.588)),
        np.array(( 4.568,  5.736,  6.905,  8.074,  10.680)),
        np.array(( 4.574,  5.796,  7.017,  8.238,  10.660)),
        np.array(( 4.568,  5.791,  7.013,  8.236,  10.681))]
statfiles = glob.glob('*.stat')
statfiles.sort()

j= 0
wantedz = [10., 6., 4., 2., 1., 0.5, 0.0] 
colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, len(wantedz))))

for i in range(len(statfiles)):
    tipsyfile = ('.').join(statfiles[i].split('.')[0:4])
    tipsy = nptipsyreader.Tipsy(tipsyfile)
    tipsy._read_param()
    tipsy._read(skipgas=True, skipdark=True, skipstar=True)
    redshift = 1./tipsy.t-1.
    if np.round(redshift, 1) in wantedz:
        stat = readstat.readstat(statfiles[i])
        nonzero = stat['starmass'] > 0.
        #nbins = 1.05*(np.log10(np.max(stat['starmass'])) - np.log10(np.min(stat['starmass'][nonzero])))

        hist, bins = np.histogram(np.log10(stat['starmass'][nonzero]*0.7), bins=nbins[j])

        print bins
        j += 1
        
        if bins is None: pdb.set_trace()
        x = (bins[1:]+bins[:-1])/2.
        normalization = (bins[1:]-bins[:-1])*(np.float(tipsy.paramfile['dKpcUnit'])/1000.*(tipsy.h/0.7))**3.
        print normalization
        print hist
        err = 1. + (hist + 0.75)**0.5
        #print "{0:.2e}".format(bins)
        thiscolor = next(colors)

        plt.errorbar(10.**x, hist/normalization, yerr = err/normalization, color=thiscolor, linewidth=linewidth, label='z~' + "{:.2}".format(redshift))
        if np.round(redshift, 1) in gonzredshift:
            ind = np.where(np.round(redshift, 1) == gonzredshift)[0]
            plt.plot(10.**gonzmasses, 10.**gonzdndm[ind], color=thiscolor, label='z~' + "{:.2}".format(gonzredshift[ind]) + ' Gonzales 2010')

        if np.round(redshift, 1) in czakredshift:
            ind = np.where(np.round(redshift, 1) == czakredshift)[0]
            plt.plot(10.**czakmasses, 10.**czakdndm[ind], color=thiscolor, label='z~' + "{:.2}".format(czakredshift[ind]) + ' Tomczak 2013')

plt.yscale('log')
plt.xscale('log')
plt.xlabel('M$_{*}$ [M$_{\odot}$]', fontsize=15)
plt.ylabel('$\phi$ [dlogM$^{-1}$ Mpc$^{-3}$]', fontsize=15)
plt.title('Cosmo6.25PLK Evolution of Stellar Mass Function')
plt.legend(ncol=2)
plt.show()
