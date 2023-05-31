#!/usr/bin/python3.10

import sys
import numpy as np
from py_spec import HMHDslab
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [10., 6]

num_args = len(sys.argv)

if(num_args<1):
    raise ValueError("Not enough input arguments!")
else:
    for n in range(1, num_args):
        HMHDslab.plot_w_vs_time(sys.argv[1])

root_name = sys.argv[1]
if(root_name[-1]=='/'):
    root_name = root_name[:-1]
plt.savefig(root_name+"_wvstime.png", dpi=400)

plt.show(block=True)
