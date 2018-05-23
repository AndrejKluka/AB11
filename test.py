import time
start = time.clock()
#----------------------------------------------------------Modules for general life
from os import path
import os
import scipy.io
import scipy
import numpy as np
import copy
from pyevtk.hl import gridToVTK
import h5py
import warnings

data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),'LHS_K_3_N_15_c_00')
data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)

print(data)
x=data['LHS']
print(x)
a=x.todense()
