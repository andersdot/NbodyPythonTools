#!/usr/bin/env python

import sys
import struct

def read_bgc2(filename):
  offset = 4
  groupoffset = 8
  particleoffset = 8

  headersize = 1024
  groupsize = 4*8 + 10*4
  particlesize = 1*8 + 6*4

  headerformat = '=Q 16q 19d'
  groupformat = '=2q 2Q 10f'
  particleformat = '=q 6f'

  print "Reading "+filename+"..."
  fd = open(filename, 'rb')
  bin_string = fd.read()
  fd.close()
  print "Finished reading file."
  bin_string = bin_string[offset:]

  # Header stuff
  header_bin = bin_string[:headersize]
  header_pad = headersize - 36*8
  header = struct.unpack(headerformat, header_bin[:-header_pad])

  # Group stuff
  ngroups = header[8]
  print 'ngroups = ', ngroups
  groupstart = headersize + groupoffset
  groupend = groupstart + ngroups*groupsize
  group_bin = bin_string[groupstart:groupend]
  group = []
  for i in range(ngroups):
    group.append(struct.unpack(groupformat, group_bin[i*groupsize:(i+1)*groupsize]))

  # Particle stuff
  particlestart = headersize + groupoffset + ngroups*groupsize + particleoffset
  particle_bin = bin_string[particlestart:]
  particle = []
  p_start = 0
  for i in range(ngroups):
    npart = group[i][2]
    particle.append([])
    for j in range(npart):
      particle[i].append(struct.unpack(particleformat, particle_bin[p_start:p_start+particlesize]))
      p_start += particlesize
    p_start += particleoffset

  return header, group, particle


def main():
  header, group, particle = read_bgc2(sys.argv[1])

  print 'Header contents:'
  for value in header:
    print value
  print

  print 'Group[0] contents:'
  for value in group[0]:
    print value
  print

  print 'Particles in group[0]:'
  for part in particle[0]:
    print part
  print

  print 'Group[1] contents:'
  for value in group[1]:
    print value
  print

  print 'Particles in group[1]:'
  for part in particle[1]:
    print part





if __name__ == '__main__':
  main()
