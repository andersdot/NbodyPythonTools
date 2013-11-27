import nptipsyreader
import glob
import numpy as np
import matplotlib.pyplot as plt

tipsyfiles = glob.glob('*.den')



totalbaryons = np.zeros(len(tipsyfiles))
coolgas      = np.zeros(len(tipsyfiles))
warmgas      = np.zeros(len(tipsyfiles))
hotgas       = np.zeros(len(tipsyfiles))
stars        = np.zeros(len(tipsyfiles))
times        = np.zeros(len(tipsyfiles))
redshift     = np.zeros(len(tipsyfiles))

for i in range(len(tipsyfiles)):
    tipsyfile = tipsyfiles[i][:-4]
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

width = 0.1
plt.bar(redshift, coolgas/totalbaryons, width, color='b', label='coolgas')
plt.bar(redshift+width, warmgas/totalbaryons, width, color='g', label='warmgas')
plt.bar(redshift+2.*width, hotgas/totalbaryons, width, color='r', label='hotgas')
plt.bar(redshift+3.*width, stars/totalbaryons, width, color='k', label='stars')

plt.legend()
plt.xlabel('time [Gyrs]')
plt.ylabel('Fraction Baryons')
plt.show()
