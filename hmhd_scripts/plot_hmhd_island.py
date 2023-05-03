#!/usr/bin/python3.10

import sys
import numpy as np
from py_spec import HMHDslab
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [7., 7]

num_args = len(sys.argv)

if(num_args<2):
    raise ValueError("Not enough input arguments!")
else:
    for n in range(1, num_args):

        HMHDslab.get_width_As_HMHD(sys.argv[n])
        plt.ylim([0, 2*np.pi])

plt.show(block=True)
