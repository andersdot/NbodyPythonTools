"""
Python scripts to create the scripts necessary
to run Rockstar on a series of ChaNGa outputs
in a directory

It expects .den files for each output and can
create .iord files for dm only simulations
which Rockstar requires to run

Subpackages
-----------
rockstar.makeiord
  Creates iord files,required by rockstar, for dm only runs

rockstar.snaps
  Creates the snaps.txt file containing all the snaps in the 
  directory that rockstar will analyze

rockstar.rockstarscript
  Creates the various scripts and config files that Rockstar
  requires to run

rockstar.prepRockstar.py
  A wrapper for all the rockstarscript modules. Takes queue, 
  server interface, rockstar main directory as inputs

rockstar.comparegenmf_simsrockstar
  Creates mass functions from the natural rockstar outputs
  and compares it to theoretical mass functions created at
  http://hmf.icrar.org/
"""

#Define rockstar main directory
rockstar_main_dir = '/astro/users/lmanders/code/Rockstar-Galaxies/'
#Compile all possible mf for various redshifts, cosmologies we use
#and put them in a data directory in this directory

from . import bgc2, makeiord, rockstarscript, readout
