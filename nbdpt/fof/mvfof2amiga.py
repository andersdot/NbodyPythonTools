import glob
import os
import shutil

def move():
    filenames = glob.glob('*.grp')
    for i in filenames:
        prefix = '.'.join(i.split('.')[0:3])
        shutil.copyfile(i, prefix + '.amiga.grp')
        
