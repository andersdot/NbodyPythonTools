import os
import re

def find(basedir='.'):

    return [f for f in os.listdir(basedir) if re.match('^(cosmo|h).*\.[\d]*$', f)]
