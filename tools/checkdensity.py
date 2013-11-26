import glob
import nptipsyreader
import numpy as np
import pdb
import matplotlib.pyplot as plt
import matplotlib as mpl

def averageden():
    gtpfiles = glob.glob('*.gtp')
    gtpfiles.sort()
    avgden = np.zeros(len(gtpfiles), dtype='float')
    medianden = np.zeros(len(gtpfiles), dtype='float')
    a = np.zeros(len(gtpfiles), dtype='float')
    tipsyfile = ('.').join(gtpfiles[0].split('.')[0:4])
    tipsy = nptipsyreader.Tipsy(tipsyfile)
    tipsy._read_param()
    plt.clf()
    colors = iter(mpl.cm.gist_rainbow(np.linspace(0, 1, 13)))#len(gtpfiles))))
    for i in range(len(gtpfiles)):
        gtp = nptipsyreader.Tipsy(gtpfiles[i])
        gtp._read()
        mass = gtp.star['mass']
        radius = gtp.star['eps']
        den = mass/(4/3.*np.pi*radius**3)
        avgden[i] = np.mean(den)
        medianden[i] = np.median(den)
        a[i] = gtp.t
        pdb.set_trace()
        if (gtp.t > 1/2.):
            highmass = mass*np.float(tipsy.paramfile['dMsolUnit']) > 1.
            plt.scatter(mass[highmass]*np.float(tipsy.paramfile['dMsolUnit']), den[highmass], color=next(colors),label='{:.2f}'.format(gtp.t))
            #plt.hist(den, color=next(colors), bins=20., histtype='step', label='{:.2f}'.format(gtp.t), log=True)
    plt.xscale('log')
    plt.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.5)
    plt.xlabel('M$_{h}$ [M$_{\odot}$]', fontsize=15)#('density [rho crit]', fontsize='large')
    plt.ylabel('rho crit', fontsize=15) #('logN', fontsize='large')
    plt.title('cosmo6.25PLK Halo Densities', fontsize='large')
    plt.show()
    plt.savefig('densities_mass.png')



#plt.plot(a, beta)
#plt.plot(a, avgden/(beta/0.25)**3.)
#return avgden, medianden
