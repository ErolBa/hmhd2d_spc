#!/usr/bin/env /home/balkovic/virtualenvs/venv-intel/bin/python

from py_spec import input_dict, HMHDslab
import sys
import subprocess
from textwrap import dedent
from glob import glob
import os
import json
from filelock import FileLock

new_tmax = int(sys.argv[1])
rsifile = sys.argv[2]
root_name = os.getcwd().split('/')[-1]

if(rsifile=="latest"):
    rsfiles = glob("rsTearing*.00000")
    rsifile = sorted(rsfiles, key=os.path.getmtime)[-1][:-6]
    print("restarting from", rsifile)

HMHDslab.set_inputfile_var('Tearing', 'tmax', new_tmax)
HMHDslab.set_inputfile_str('Tearing', 'rsifile', rsifile)

files=glob(rsifile+".*")

# get the inputs dict
inputs = input_dict()
json_fname = glob("inputs*")[0]
with open(json_fname, 'r') as fh:
    temp_dict = json.load(fh)
inputs.update(temp_dict)

# lock file so only one code changes files and runs make at the time
lock = FileLock("/home/balkovic/codes/hmhd2d_spc/build/lock.lock")
with lock:

    HMHDslab.set_hmhd_grid(inputs.mesh_pts[1], inputs.mesh_pts[0], inputs.num_cpus[1], inputs.num_cpus[0], "~/codes/hmhd2d_spc")    
    subprocess.run(f"cd ~/codes/hmhd2d_spc/build; make", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)     
    
    with open("slurmjobrst.run", "w") as f:
        print(dedent(f"""\
            #!/bin/sh
            #SBATCH --nodes 1
            #SBATCH --ntasks {len(files)}
            #SBATCH --cpus-per-task 1
            #SBATCH --mem 8G
            #SBATCH --time 01:00:00
            #SBATCH --job-name="hmhd2d_rst {root_name}"
            export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
            srun /home/balkovic/codes/hmhd2d_spc/build/hmhd2d Tearing"""), 
            file=f)

proc = subprocess.Popen(f"sbatch -W slurmjobrst.run", shell=True)

# wait for the slurm job to finish
proc.wait()
