#!/usr/bin/env /home/balkovic/virtualenvs/venv-intel/bin/python

from glob import glob
from re import match
from subprocess import run

run(f"squeue --user=balkovic -o \" %j\" -h > running_jobs.txt", shell=True)
with open("running_jobs.txt", "r") as f:
    print("RUNNING JOBS:")
    for l in f:

        # latest_dir = max(glob('*'), key=os.path.getmtime)        
        latest_dir = l.split(" ")[2][:-1]

        len_files = len(glob(latest_dir+"/data*.hdf"))

        with open(latest_dir+"/iTearing", "r") as f:
            for line in f:
                tpltxint = match('^tpltxint.*', line)
                if(tpltxint is not None):
                    tpltxint = tpltxint.group(0)
                    break

        t_interval = float(tpltxint.split('=')[1].split('!')[0].replace("D", "E"))
        time_completed = (len_files-1) * t_interval

        with open(latest_dir+"/iTearing", "r") as f:
            for line in f:
                time_total = match('^tmax.*', line)
                if(time_total is not None):
                    time_total = time_total.group(0)
                    break

        time_total = float(time_total.split('=')[1].split('!')[0].replace("D", "E"))

        print(f"\tFor run '{latest_dir}' completed {time_completed:.4e} / {time_total:.4e} alfven times")

run("rm running_jobs.txt", shell=True)