"""
A set of tools to help analyze the simulations

Modules
--------

checksimlog
  parses through the log file while a simulation is running
  and plots time/step v time [Gyrs] where each restart is a 
  different color

cosmography
  some simple (but needing update) programs to calculate
  cosmological time, distance for a given redshift

fixdecomp
  will create a tipsy array file to read into tipsy/salsa to mark 
  particles of each classification a different color

nptipsyreader
   The reason I started this project. Pynbody originally wasn't compiling
   on some of the supercomputers so I wrote this code as a way to read in 
   sims to get units and such quickly. It reads large files in wicked fast
   and reads in the param file, calculating some units for you

readstat
  Expects a certain structure of the stat file, but will read it in. I need
  to add additional options for the strutcture of the stat file depending
  on which halo finder created it. Also Rockstar's natural outputs

starlog
  Reads in a starlog file, and converts units to Msol, Gyr

Subpackages
----------
analysis
  various plotting routines

rockstar
  scripts to build scripts to run rockstar

fof
  scripts associated with running a fof algorithm on a simulation

"""

from . import nptipsyreader, readstat, starlog, statfile, checksimlog, cosmography

from . import analysis, rockstar, fof
