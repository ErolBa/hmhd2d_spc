#!/usr/bin/env python3

import sys
import numpy as np
from py_spec import HMHDslab
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [7., 7]

num_args = len(sys.argv)

if(num_args<3):
    raise ValueError("Not enough input arguments!")
else:
    for n in range(2, num_args):
        HMHDslab.plot_scalar(sys.argv[n], sys.argv[1], ylims=[-np.pi, np.pi])
        plt.tight_layout()

plt.show(block=True)
