import glob

files = glob.glob('*.iord')
files.sort()
snaps = []
f = open('snaps.txt', 'w')

for i in files: f.write(i.split('.')[-2] + '\n')

f.close()
