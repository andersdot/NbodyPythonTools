"""
The basic functions to create usable rockstar scripts to run rockstar on all tipsy
outputs in a directory. 

snaps()
  create a snaps.txt file containing all the outputs in the directory
  one number per line

cfg(ncorespernode=32, nnodes=1, ServerInterface='ipogif0', massdef='200c', massdef2=None)
  create a rockstar.submit.cfg file to submit to the rockstar program
  nocorespernode: the number of cores on each node of the computer
  ServerInterface: 'ib0' for most computers except bluewaters which uses 'ipogif0'
  massdef: how to define the virial mass of a halo 
  massdef2: if you want another mass defined

mainsubmissionscript(walltime = '24:00:00', email = 'l.sonofanders@gmail.com', machine='stampede', nnodes=1, ncorespernode=32, queue='largemem', rockstardir='/home1/02575/lmanders/code/Rockstar-Galaxies/')
  create the submission script to run rockstar. 
  walltime: total amount of wallclock time alloted to the analysis
  email: the email you want messages from the queue sent to
  machine: the machine you will run the process on; 
           right now either 'stampede', 'pleiades' or 'bluewaters'
  nnodes: the number of nodes to run on
  ncorespernode: the number of cores per node on the machine
  queue: the queue to submit to
  rockstardir: the directory your rockstar copy lies in

postsubmissionscript(email = 'l.sonofanders@gmail.com', machine = 'stampede', queue = 'largemem'\
, rockstardir='/home1/02575/lmanders/code/Rockstar-Galaxies', ncorespernode=32, walltime='24:00:00')
  create the submission script to do post processing on rockstar and turn its natural outputs 
  into .grp and .stat files
  email: the email you want messages from the queue sent to
  machine: the machine you will run the process on; 
           right now either 'stampede', 'pleiades', or 'bluewaters' 
  queue: the queue to submit to 
  rockstardir: the directory your rockstar copy lies in 
  ncorespernode: the number of cores per node on the machine 
  walltime: total amount of wallclock time alloted to the analysis
"""


import glob
import numpy as np
import struct
from .. import nptipsyreader
import os
from .. import findtipsy

def snaps():
    files = glob.glob('*.iord')
    files.sort()
    snaps = []
    f = open('snaps.txt', 'w')
    for i in files: f.write(i.split('.')[-2] + '\n')
    f.close()

def cfg(ncorespernode=32, nnodes=1, ServerInterface='ipogif0', massdef='200c', massdef2=None):
    filename = ('.').join(glob.glob('*.iord')[0].split('.')[:-2])
    tipsyfile = ('.').join(glob.glob('*.iord')[0].split('.')[:-1])
    tipsy = nptipsyreader.Tipsy(tipsyfile)
    tipsy._read_param()
    dmsol = np.float(tipsy.paramfile['dMsolUnit'])
    dkpc = np.float(tipsy.paramfile['dKpcUnit'])
    f = open(tipsyfile)
    (t, nbodies, ndim, nsph, ndark, nstar) = struct.unpack('>diiiii', f.read(28))
    f.seek(struct.calcsize('>diiiiii')+struct.calcsize('12f4')*nsph)
    (mass, x, y, z, vx, vy, vz, eps, phi) = struct.unpack('>9f4', f.read(36))
    darkmass = mass*dmsol*tipsy.h
    #f.seek(struct.calcsize('>diiiiii')+struct.calcsize('12f4')*nsph+struct.calcsize('9f4')*ndark)
    #(mass,x,y,z,vx,vy,vz,metals,tform,eps,phi) = struct.unpack('>11f4', f.read(44))
    force = eps*dkpc*tipsy.h/1e3
    f.close()
    
    f = open('rockstar.submit.cfg', 'w')
    f.write('OVERLAP_LENGTH=0.5 \n')
    f.write('PARALLEL_IO=1 \n')
    f.write('NUM_WRITERS=' + str(ncorespernode*nnodes)+' \n')
    f.write('NUM_BLOCKS=1 \n')
    f.write('FILENAME=' + filename +'.<snap> \n')
    f.write('#NUM_SNAPS=1 \n')
    f.write('SNAPSHOT_NAMES=snaps.txt  #One snapshot name per line \n')
    f.write('FORK_READERS_FROM_WRITERS=1 \n')
    f.write('FORK_PROCESSORS_PER_MACHINE=' + str(ncorespernode) + '\n')
    f.write('FILE_FORMAT = "TIPSY" # or "ART" or "ASCII" \n')
    f.write('PARTICLE_MASS ='+ str(darkmass) +'  #Msol/h \n')
    f.write('FORCE_RES = ' + str(force) + '  #Mpc/h \n')                                            
    f.write('# You should specify cosmology parameters only for ASCII formats \n# For GADGET2 and ART, these parameters will be replaced with values from the \n# particle data file \n')
    f.write('h0 =' + str(tipsy.h) + '\n')
    f.write('Ol =' + tipsy.paramfile['dLambda'] + '\n')
    f.write('Om =' + tipsy.paramfile['dOmega0'] + '\n')
    f.write('# For GADGET2, you may need to specify conversion parameters. \n# Rockstar"s internal units are Mpc/h (lengths) and Msun/h (masses) \n')                              
    f.write('TIPSY_MASS_CONVERSION   = ' + str(dmsol*tipsy.h) + '\n')
    f.write('TIPSY_LENGTH_CONVERSION = '  + str(dkpc*tipsy.h/1e3) + '\n')
    f.write('TIPSY_VELOCITY_CONVERSION = ' + str(tipsy.velunit) + '\n')
    f.write('BOX_SIZE = ' + str(dkpc*tipsy.h/1e3) + '\n')
    f.write('#For ascii files, the file format is assumed to be: \n# X Y Z VX VY VZ ID \n')
    f.write('#For non-periodic boundary conditions only: \n#PERIODIC=0 \n')                    
    f.write('FULL_PARTICLE_CHUNKS = ' + str(ncorespernode*nnodes)  + '  # Should be the same as NUM_WRITERS above to save all particles \n')
    f.write('PARALLEL_IO_SERVER_INTERFACE = "' + ServerInterface + '"\n')
    f.write('MASS_DEFINITION  = "' + massdef + '"\n')
    if massdef2: f.write('MASS_DEFINITION2 = "' + massdef2+ '"\n')

def mainsubmissionscript(walltime = '24:00:00', email = 'l.sonofanders@gmail.com', machine='stampede', nnodes=1, ncorespernode=32, queue='largemem', rockstardir='/home1/02575/lmanders/code/Rockstar-Galaxies/'):
    sbatchname = findtipsy.find()[0].split('.')[0]
    if (machine == 'pleiades') or (machine == 'bluewaters'):
        filename = 'rockstar.qsub'
        f = open(filename, 'w')

        f.write('#!/bin/bash \n')
        f.write('#PBS -N ' + sbatchname + ' \n')
        if machine == 'pleiades': f.write('#PBS -lselect=' + str(nnodes) + ':ncpus='+str(ncorespernode)+':mpiprocs=1'+ '\n')
        if machine == 'bluewaters':f.write('#PBS -lnodes=' + str(nnodes) + ':ppn='+str(ncorespernode)+'\n')
        f.write('#PBS -m be \n')
        f.write('#PBS -M ' + email +' \n')
        f.write('#PBS -q ' + queue + ' \n')
        f.write('#PBS -l walltime=' + walltime + ' \n')
        
        f.write('cd $PBS_O_WORKDIR \n')
        f.write('rm auto-rockstar.cfg \n')
        f.write('ulimit -c unlimited \n')

        if machine=='bluewaters': prefix = 'aprun -n ' + str(nnodes) + ' -N 1 '
        else: prefix = ''
    
        f.write(rockstardir + 'rockstar-galaxies -c rockstar.submit.cfg & \n')
        f.write("perl -e 'sleep 1 while (!(-e "+'"' + "auto-rockstar.cfg"+'"' + "))' \n")
        f.write(prefix + rockstardir + 'rockstar-galaxies -c auto-rockstar.cfg \n')

    if machine == 'stampede':
        filename = 'rockstar.sbatch'
        f = open(filename, 'w')

        f.write('#!/bin/bash \n')
        f.write('#SBATCH -J ' + sbatchname + ' \n')
        f.write('#SBATCH -n ' + str(nnodes) + ' \n')
        f.write('#SBATCH -N ' + str(nnodes) + ' \n')
        f.write('#SBATCH -o ' + sbatchname + '.rockstar.o%j \n')
        f.write('#SBATCH -p ' + queue + ' \n')
        f.write('#SBATCH -t 24:00:00 \n')
        f.write('#SBATCH --mail-user=l.sonofanders@gmail.com \n')
        f.write('#SBATCH --mail-type=ALL \n')
        f.write('#SBATCH -A TG-MCA94P018 \n')
        f.write('## -n: total tasks \n')
        f.write('## -N: nodes \n')
        
        f.write('cd $SLURM_SUBMIT_DIR \n')
        f.write('rm auto-rockstar.cfg \n')
        f.write('ulimit -c unlimited \n')

        f.write(rockstardir + 'rockstar-galaxies -c rockstar.submit.cfg & \n')
        f.write("perl -e 'sleep 1 while (!(-e "+'"' + "auto-rockstar.cfg"+'"' + "))' \n")
        f.write('ibrun ' + rockstardir + 'rockstar-galaxies -c auto-rockstar.cfg \n')



def postsubmissionscript(email = 'l.sonofanders@gmail.com', machine = 'stampede', queue = 'largemem', rockstardir='/home1/02575/lmanders/code/Rockstar-Galaxies', ncorespernode=32, walltime='24:00:00', nnodes=1):
    tipsyfile = findtipsy.find()[0]
    iordfilepre = ('.').join(tipsyfile.split('.')[:-1])
    
    
    sim = nptipsyreader.Tipsy(tipsyfile)
    sim._read_param()
    boxsize = str(float(sim.paramfile['dKpcUnit'])/1000.*sim.h)
    
    outfiles = glob.glob('out_*.list')
    
    snapsf = open('snaps.txt')
    snaps = [line.split('\n')[0] for line in snapsf]
    snapsf.close()
    
    genstatexecline = []
    
    if machine == 'stampede':
        f = open('rockstar.post.sbatch', 'w')

        f.write('#!/bin/bash \n')
        f.write('#SBATCH -J' + iordfilepre +' \n')
        f.write('#SBATCH -n 1 \n')
        f.write('#SBATCH -N 1 \n')
        f.write('#SBATCH -o ' + iordfilepre + '.rockstar.o%j \n')
        f.write('#SBATCH -p ' + queue + ' \n')
        f.write('#SBATCH -t ' + walltime + ' \n')
        f.write('#SBATCH --mail-suers=' + email + '\n')
        f.write('#SBATCH --mail-type=ALL \n')
        f.write('#SBATCH -A TG-MCA94P018 \n')
        f.write('cd $SLURM_SUBMIT_DIR \n')
        f.write('rm auto-rockstar.cfg \n')
        f.write('ulimit -c unlimited \n')
    if (machine == 'pleiades') or (machine == 'bluewaters'):
        f = open('rockstar.post.qsub', 'w')

        f.write('#!/bin/bash \n')
        f.write('#PBS -N cosmo6.rockstar \n')
        if machine == 'pleiades':f.write('#PBS -lselect=1:ncpus='+str(ncorespernode)+':mpiprocs=1 \n')
        if machine == 'bluewaters':f.write('#PBS -lnodes=1:ppn='+str(ncorespernode)+' \n')
        f.write('#PBS -q ' + queue + ' \n')
        f.write('#PBS -l walltime=24:00:00 \n')
        f.write('#PBS -m be \n')
        f.write('#PBS -M ' + email + '\n')
        f.write('cd $PBS_O_WORKDIR \n')
        f.write('ulimit -c unlimited \n')

    

    for i in range(len(snaps)):    
        parentexecline = (rockstardir + 'util/find_parents out_'+ str(i  ) + '.list ' + boxsize + ' > out_'+ snaps[i] + '.parents \n' )
        f.write(parentexecline)
        genstatexecline.append(rockstardir + 'examples/gen_grp_stats out_' + snaps[i] + '.parents ' + iordfilepre + '.' +snaps[i] + '.iord  halos_' + snaps[i] + '.*.particles \n')
        
        
    for i in range(len(snaps)):
        #f.write("perl -e 'sleep 1 while (!(-e " + '"' + "out_" + snaps[i] + ".parents"+'"'+"))'" + "\n")
        f.write(genstatexecline[i])
    for i in range(len(snaps)):
        f.write("perl -e 'sleep 1 while (!(-e " + '"' + "out_" + str(i) + ".grp"+'"'+"))'" + "\n")
        f.write("mv out_" + str(i) + ".grp "  + iordfilepre+'.'+snaps[i]+".rockstar.grp \n")
        f.write("perl -e 'sleep 1 while (!(-e " + '"' + "out_" + str(i) + ".stat"+'"'+"))'" + "\n")
        f.write("mv out_" + str(i) + ".stat " + iordfilepre+'.'+snaps[i]+".rockstar.stat \n")
    f.close()
    #    os.system(parentexecline)
    #    os.system(genstatexecline)
