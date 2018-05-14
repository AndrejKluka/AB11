from __future__ import print_function
from os import path
import scipy.io
import numpy as np
import time
from math import *
import pprint
import pylab
import matplotlib
import matplotlib.pyplot as plt

start = time.time()
# reading raw dataset and putting them into u,v,w arrays
##vspace_dir=path.join(path.dirname(_file_),'vspace.npy') 
data=scipy.io.loadmat('raw_data_1', mdict=None, appendmat=True)
#listx=np.loadtxt(path.join(txt_dir,'xlist.txt')) 
u=data['u']
w=data['w']
v=data['v']

def mean(values):

      length = len(values)
      tot_sum = 0

      for i in range(length):
            tot_sum += values[i]
      tot_sum = sum(values)
      tot_sum2 = sum(tot_sum)
      tot_sum3 = sum(tot_sum2)
      average = tot_sum3/length**3
      return average


m = mean(u)
n = mean (v)
q = mean (w)
print("Mean of u:", m)
print("Mean of v:", n)
print("Mean of w:", q)

#------------------------------------------------------------#

print("Variance of u:", np.var(u[0][0]))
print("Variance of v:", np.var(v[0][0]))
print("Variance of w:", np.var(w[0][0]))


end = time.time()
print("Time taken to compute:", end - start)
