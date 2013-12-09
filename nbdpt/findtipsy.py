import os
import re

def find():

    return [f for f in os.listdir('.') if re.match('^[\d]*$', f.split('.')[-1])]
