import numpy as np

def read(filename):
    f = open(filename)
    names = f.readline().strip('\n#').split(' ')
    dtype = np.concatenate((['int', 'int'], (36*'float ').split(' ')[:-1]))
    f.close()
    stat = np.genfromtxt(filename, dtype=dtype, names=names, skiprows=1)
    return stat
