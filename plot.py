from __future__ import division
import mdtraj as md
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import numpy as np
# from numpy import zeros, sqrt, square, log, where, pi, mean, average, arange, sin, cos
from math import sqrt, log, pi, sin, cos
from textwrap import wrap
from cycler import cycler
from glob import glob
import sys
from itertools import groupby
from decimal import Decimal
from mpmath import mp
import random

data = glob('data/*')
print data

for i in data:
    print i
    t = np.load(i)
    print t
    print t.shape
    print