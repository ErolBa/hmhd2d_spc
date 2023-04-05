#!/usr/bin/env python3

import sys
import numpy as np
from py_spec import HMHDslab
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [7., 7]

HMHDslab.plot_scalar(sys.argv[1], sys.argv[2], xlims=[-2.5, 2.5], ylims=[-np.pi, np.pi])

plt.tight_layout()
plt.show(block=True)
