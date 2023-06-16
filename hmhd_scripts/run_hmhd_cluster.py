#!/usr/bin/env /home/balkovic/virtualenvs/venv-intel/bin/python

from py_spec import input_dict, HMHDslab
import json
import sys
import subprocess
from filelock import FileLock
from textwrap import dedent

# get the inputs dict
inputs = input_dict()
json_fname = sys.argv[1]
with open(json_fname, 'r') as fh:
    temp_dict = json.load(fh)
inputs.update(temp_dict)

# run HMHD
inputs.config = HMHDslab.gen_profiles_from_psi(inputs.psi_profile)

subprocess.run(f"cp ~/codes/hmhd2d_spc/iTearing .", shell=True)
HMHDslab.set_inputfile_var('Tearing', "tmax", inputs.tmax)
HMHDslab.set_inputfile_var('Tearing', "tpltxint", inputs.tpltxint)
HMHDslab.set_inputfile_var('Tearing', "xl", inputs.xl)
HMHDslab.set_inputfile_var('Tearing', "dt", inputs.dt) # d is 2e-3
HMHDslab.set_inputfile_var('Tearing', "eta", inputs.eta)
HMHDslab.set_inputfile_var('Tearing', "visc", inputs.visc)
HMHDslab.set_inputfile_var('Tearing', "qpb", inputs.qpb)
HMHDslab.set_inputfile_var('Tearing', "qpc", inputs.qpc)
HMHDslab.set_inputfile_var('Tearing', "bs_curr_const", inputs.bs_curr_const)

# lock file so only one code changes files and runs make at the time
lock = FileLock("/home/balkovic/codes/hmhd2d_spc/build/lock.lock")
with lock:
    HMHDslab.set_hmhd_dens(inputs.dens_profile, "~/codes/hmhd2d_spc")
    HMHDslab.set_hmhd_pres(inputs.pres_profile, "~/codes/hmhd2d_spc")
    HMHDslab.set_hmhd_profiles(inputs.config, "~/codes/hmhd2d_spc")
    HMHDslab.set_hmhd_heatcond(inputs.heatcond_perp, inputs.heatcond_para, inputs.heatcond_flag, "~/codes/hmhd2d_spc")
    
    HMHDslab.set_hmhd_grid(inputs.mesh_pts[1], inputs.mesh_pts[0], inputs.num_cpus[1], inputs.num_cpus[0], "~/codes/hmhd2d_spc")
    
    subprocess.run(f"cd ~/codes/hmhd2d_spc/build; make", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)     
    subprocess.run(f"cp ~/codes/hmhd2d_spc/build/globals.f90 .; cp ~/codes/hmhd2d_spc/build/prob.f90 .; cp ~/codes/hmhd2d_spc/build/equation.f90 .", shell=True)

    with open("slurmjob.run", "w") as f:
        print(dedent(f"""\
            #!/bin/sh
            #SBATCH --nodes 1
            #SBATCH --ntasks {inputs.num_cpus[0]*inputs.num_cpus[1]}
            #SBATCH --cpus-per-task 1
            #SBATCH --mem 8G
            #SBATCH --time 01:00:00
            #SBATCH --job-name="hmhd2d {inputs.root_fname}"
            export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
            srun /home/balkovic/codes/hmhd2d_spc/build/hmhd2d Tearing"""), 
            file=f)
    
    proc = subprocess.Popen(f"sbatch -W slurmjob.run", shell=True)

# wait for the slurm job to finish
proc.wait()
