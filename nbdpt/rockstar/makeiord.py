from .. import nptipsyreader
import glob
import shutil
import time

def make():
    tipsyfiles = glob.glob('*.den')
    tipsyfiles.sort()
    sim = nptipsyreader.Tipsy(tipsyfiles[0][:-4])
    (t, nbodies, ndim, nsph, ndark, nstar) = sim.unpack_header()
    print t, nbodies, ndim, nsph, ndark, nstar
    f = open(tipsyfiles[0][:-4] + '.iord', 'w')
    f.write(str(nbodies) + '\n')
    for i in range(nbodies): f.write(str(i) + '\n')
    f.close()
    for i in tipsyfiles[1:]: 
        shutil.copy(tipsyfiles[0][:-4] + '.iord', i[:-4] + '.iord')
        time.sleep(2)
