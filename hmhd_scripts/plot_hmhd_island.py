#!/usr/bin/python3.10

import sys
import numpy as np
from py_spec import HMHDslab
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [7., 7]

HMHDslab.get_width_As_HMHD(sys.argv[1])

#plt.ylim([np.pi-0.8, np.pi+0.8])
plt.ylim([0, 2*np.pi])

plt.show(block=True)
