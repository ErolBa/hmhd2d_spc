#!/usr/bin/env python3

import sys
import matplotlib
matplotlib.use("Qt5Agg")
from py_spec import HMHDslab

HMHDslab.plot_scalar_3d(sys.argv)