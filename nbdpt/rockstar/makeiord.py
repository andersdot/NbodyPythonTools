from .. import nptipsyreader
import glob
import shutil
import time
import os
import re
from .. import findtipsy

def make():
    tipsyfiles = findtipsy.find()
    tipsyfiles.sort()
    sim = nptipsyreader.Tipsy(tipsyfiles[0])
    (t, nbodies, ndim, nsph, ndark, nstar) = sim.unpack_header()
    print t, nbodies, ndim, nsph, ndark, nstar
    f = open(tipsyfiles[0] + '.iord', 'w')
    f.write(str(nbodies) + '\n')
    for i in range(nbodies): f.write(str(i) + '\n')
    f.close()
    for i in tipsyfiles[1:]: 
        shutil.copy(tipsyfiles[0] + '.iord', i + '.iord')
        time.sleep(2)
