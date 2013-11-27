import foflinkinglength
from .. import nptipsyreader
import glob
import numpy as np
import pdb

def writescript():
    betacoeff = foflinkinglength.getbeta()
    
    denfiles = glob.glob('*.den')
    denfiles.sort()
    
    alreadymade_fof = glob.glob('*.fof.*.grp')
    alreadymade_fof.sort()
    
    found       = [('.').join(filename.split('.')[:-1]) for filename in denfiles]
    alreadymade = [('.').join(filename.split('.')[:-4]) for filename in alreadymade_fof]
    
    wanted = list(set(alreadymade) ^ set(found))
    wanted.sort()
    
    fofcode = '/u/sciteam/lmanders/code/fof/fof'
    
    fofscript = open('fof.qsub', 'w')
    
    hours = len(wanted)*2
    if hours > 24: hours = 24
    
    fofscript.write('#!/bin/bash --login \n')
    fofscript.write('#PBS -lwalltime='+str(hours)+':00:00 \n')
    fofscript.write('#PBS -lnodes=1:ppn=32 \n')
    fofscript.write('#PBS -m be \n')
    fofscript.write('#PBS -M l.sonofanders@gmail.com \n')
    fofscript.write('#PBS -N fof \n')
    fofscript.write('module load PrgEnv-cray \n')
    fofscript.write('module swap PrgEnv-cray PrgEnv-gnu \n')
    fofscript.write('export OMP_NUM_THREADS=32 \n')
    fofscript.write('cd $PBS_O_WORKDIR \n')
    
    for i in range(len(wanted)):
        tipsyfile = wanted[i]
        tipsy = nptipsyreader.Tipsy(tipsyfile)
        tipsy._read(skipgas=True, skipstar=True)
        beta = betacoeff[2] + betacoeff[1]*tipsy.t + betacoeff[0]*tipsy.t**2.
        nparticles = (tipsy.ndark)**(1/3.)*2.
        linkinglength = beta/nparticles
        line = 'aprun -n 1 -d 32 '+fofcode+ ' -e '+'{:.5e}'.format(linkinglength)+ ' -m '+'{:n}'.format(nparticles)+' -v -o '+tipsyfile+'.fof.' + '{:.3f}'.format(beta) + ' -p 1 -std < '+tipsyfile +'\n'
        fofscript.write(line)
        
        fofscript.close()
        
#beta = betacoeff[2] + betacoeff[1]*a + betacoeff[2]*a**2.
if __name__ == '__main__':
    writescript()
