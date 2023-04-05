# tar all restart files in the directory

import glob
import os

# get restart file sets

sets = glob.glob('*.00000')

for name in sets:
    fhead = name[:-5]
    fname = fhead + '*'
    tarfile = 'tar_' + fhead + 'tar'
    cmd = 'tar -cvf ' + tarfile + ' ' + fname
    print(cmd)
    os.system(cmd)
    cmd = 'rm ' + fname
    os.system(cmd)
