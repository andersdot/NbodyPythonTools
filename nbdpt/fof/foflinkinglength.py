import glob
from .. import nptipsyreader
import numpy as np
import pdb
import matplotlib.pyplot as plt

def getbeta():
    gtpfiles = glob.glob('/scratch/sciteam/lmanders/cosmo6.25PLK.192g/FofOct30/*.gtp')
    gtpfiles.sort()
    avgden = np.zeros(len(gtpfiles), dtype='float')
    medianden = np.zeros(len(gtpfiles), dtype='float')
    a = np.zeros(len(gtpfiles), dtype='float')
    
    
    for i in range(len(gtpfiles)):
        gtp = nptipsyreader.Tipsy(gtpfiles[i])
        gtp._read()
        mass = gtp.star['mass']
        radius = gtp.star['eps']
        den = mass/(4/3.*np.pi*radius**3)
        avgden[i] = np.mean(den)
        medianden[i] = np.median(den)
        a[i] = gtp.t
            
    slope = (avgden[-1] - avgden[0])/(a[-1] - a[0])
    b     = avgden[0] - slope*a[0]
    
    crit_factor = 200.
    beta = 0.125*(avgden/crit_factor)**(1/3.)
    
    coeff = np.polyfit(a, beta, 2)
    betafit = coeff[2] + coeff[1]*a + coeff[0]*a**2.
    return coeff


#plt.plot(a, beta)
#plt.plot(a, avgden/(beta/0.25)**3.)
#return avgden, medianden
