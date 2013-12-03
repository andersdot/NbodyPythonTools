from distutils.core import setup

DESCRIPTION = 'Tools to help analyze large sph nbody simulations'
LONG_DESCRIPTION = open('README.md').read()
NAME = 'NbodyPythonTools'
VERSION = '1.0'
AUTHOR = 'Lauren Anderson'
AUTHOR_EMAIL = 'l.sonofanders@gmail.com'
MAINTAINER = 'Lauren Anderson'
MAINTAINER_EMAIL = 'l.sonofanders@gmail.com'
URL = 'http://salomanders.github.com/NbodyPythonTools'
DOWNLOAD_URL = 'http://github.com/salomanders/NbodyPythonTools'
LICENSE = 'BSD'



setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR, 
      author_email= AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      package_dir={'nbdpt/':''},
      packages=['nbdpt', 
                'nbdpt/fof', 
                'nbdpt/analysis', 
                'nbdpt/rockstar',
                'nbdpt/amiga'],
      package_data={'nbdpt/rockstar':['data/mf_planck13.dat',
                                      'data/mf_wmap3.dat',
                                      'data/mf_wmap7.dat']},
      classifiers = ["Development Status :: 3 - Alpha",
                     "Intended Audience :: Developers",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                     "Programming Language :: Python :: 2",
                     "Topic :: Scientific/Engineering :: Astronomy"])
