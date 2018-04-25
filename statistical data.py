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


#------------------------------------------------------------#
n = 0.85
array1 = [u[0][0][0],v[0][0][0],w[0][0][0]]
correlationu = []
correlationv = []
correlationw = []
corrsum = 0

##for j in range(0,191):
##    for p in range(0,191):
##        for r in range(0,191):
##            corrmatrix = correlation(array1,[u[j][p][r],v[j][p][r],w[j][p][r]])
##            correlation.append(corrmatrix[1][0])
##            corrsum = corrsum + corrmatrix[1][0]
##            corrmean = corrsum/191**3


for i in range(0,192):
      for j in range(0,192):
            corrmatrixu = np.correlate(array1, [u[j][i][i],v[j][i][i],w[j][i][i]])
            corrmatrixu = list(corrmatrixu)
            correlationu.append(corrmatrixu)
            correlationarrayu = np.array(correlationu)
            uterm1 = correlationu[0]
            finalcorrelationu = correlationarrayu/uterm1

      for p in range(0,192):
            corrmatrixv = np.correlate(array1, [u[i][p][i],v[i][p][i],w[i][p][i]])
            corrmatrixv = list(corrmatrixv)
            correlationv.append(corrmatrixv)
            correlationarrayv = np.array(correlationv)
            vterm1 = correlationv[0]
            finalcorrelationv = correlationarrayv/vterm1
      for r in range(0,192):
            corrmatrixw = np.correlate(array1, [u[i][i][r],v[i][i][r],w[i][i][r]])
            corrmatrixw = list(corrmatrixw)
            correlationw.append(corrmatrixw)
            correlationarrayw = np.array(correlationw)
            wterm1 = correlationw[0]
            finalcorrelationw = correlationarrayw/wterm1




meancorrelation = (finalcorrelationu + finalcorrelationv + finalcorrelationw)/3

print(finalcorrelationu)
print("Done with correlation of u")

print(finalcorrelationv)
print("Done with correlation of v")

print(correlationw)
print("Done with correlation of w")

print (meancorrelation)



plt.plot(meancorrelation)
plt.show()
print("Done with plotting")




                  

##def corrv():
##      countv=0
##      for k in range(0,191):
##            correlation = np.corrcoef(v[k][k],v[k+1][k+1])
##            if correlation[1][0] < n:
##                  ##countv += 1
##                  print(k)
##      return (countv)
##
##def corrw():
##      countw=0
##      for l in range(0,191):
##            correlation = np.corrcoef(w[l][l],w[l+1][l+1])
##            if correlation[1][0] < n:
##                  countw = countw + 1
##                  print(l)
##      return (countw)
##
##
##
##print(corru())
##print(corrv())
##print(corrw())
##
##def countu():
##      count_u=0
##      for a in range(0,191):
##            correlation = np.corrcoef(u[a][a],u[a+1][a+1])
##            if correlation[1][0] < n:
##                  count_u = count_u + 1
##      return (count_u)
##
##def countv():
##      count_v=0
##      for b in range(0,191):
##            correlation = np.corrcoef(v[b][b],v[b+1][b+1])
##            if correlation[1][0] < n:
##                  count_v = count_v + 1
##      return (count_v)
##
##def countw():
##      count_w=0
##      for c in range(0,191):
##            correlation = np.corrcoef(w[c][c],w[c+1][c+1])
##            if correlation[1][0] < n:
##                  count_w = count_w + 1
##      return (count_w)
##
##count_u = countu()
##count_v = countv()
##count_w = countw()
##totcount = count_u + count_v + count_w
##percent = (totcount*100.)/(191.*3)
##print("Percentage uncorrelated", percent)

##h = u[0:191][0][0]
##print (np.correlate(h,h, "full"))
##plt.plot(h)
##plt.show()



##print(np.corrcoef(u[191][191][191],u[189][189][189]))

end = time.time()
print("Time taken to compute:", end - start)
