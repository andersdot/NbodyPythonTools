import glob
from .. import nptipsyreader
import numpy as np
import matplotlib.pyplot as plt
from .. import findtipsy

def plot():
    tipsyfiles = findtipsy.find()
    tipsyfiles.sort()
    totalbaryons = np.zeros(len(tipsyfiles))
    coolgas      = np.zeros(len(tipsyfiles))
    warmgas      = np.zeros(len(tipsyfiles))
    hotgas       = np.zeros(len(tipsyfiles))
    stars        = np.zeros(len(tipsyfiles))
    times        = np.zeros(len(tipsyfiles))
    redshift     = np.zeros(len(tipsyfiles))
    
    for i in range(len(tipsyfiles)):
        tipsyfile = tipsyfiles[i]
        tipsy = nptipsyreader.Tipsy(tipsyfile)
        tipsy._read_param()
        tipsy._read()
        dmsol = np.float(tipsy.paramfile['dMsolUnit'])    
        times[i] = np.max(tipsy.star['tform'])*tipsy.timeunit
        redshift[i] = 1./tipsy.t -1.
        
        totalbaryons[i] = (np.sum(tipsy.gas['mass']) + np.sum(tipsy.star['mass'])) * dmsol
        coolgas[i]      = np.sum(tipsy.gas['mass'][tipsy.gas['temp'] <= 1e4]) * dmsol
        warmgas[i]      = np.sum(tipsy.gas['mass'][(tipsy.gas['temp'] <= 1e5) & (tipsy.gas['temp'] > 1e4)]) * dmsol
        hotgas[i]       = np.sum(tipsy.gas['mass'][(tipsy.gas['temp'] >1e5)]) * dmsol
        stars[i]        = np.sum(tipsy.star['mass'])*dmsol
        
    width = np.linspace(0.1, 0.01, len(redshift))
    plt.bar(redshift+1, coolgas/totalbaryons, width, color='b', label='coolgas')
    plt.bar(redshift+width+1, warmgas/totalbaryons, width, color='g', label='warmgas')
    plt.bar(redshift+2.*width+1, hotgas/totalbaryons, width, color='r', label='hotgas')
    plt.bar(redshift+3.*width+1, stars/totalbaryons, width, color='k', label='stars')
    
    plt.legend()
    plt.xlabel('log(z+1)')
    plt.xscale('log')
    plt.ylabel('Fraction Baryons')
    plt.show()

if __name__=='__main__': plot()
