# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 10:50:00 2018

@author: Lisa
"""
import time
start = time.clock()
#   Modules used for plotting 
#from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
##   import matplotlib.pyplot
#import plotly.plotly as py
#import plotly
#from plotly.graph_objs import *
#import plotly.figure_factory

from os import path
import scipy.io
import scipy
import numpy as np
from skimage import measure
import math

stop= time.clock()
print ('\n',int((stop-start)*1000)/1000.,'sec -- imported modules')

#   General setup for program run
to_load=False          # if true will load already the last calculated Q or lambda dataset
to_plotly=False        # if true will send the plot to plotly website
to_matplot=False        # if true will use matplotlib to plot
n_elements=60         # number of elements on each side of cube calculated
to_calc_Q=False         # if true will calc Q on cube with n_elements
to_calc_Lambda2=True   # if true will calc lambda2 on cube with n_elements
q_threshold=0.16 
      
data_num=1              # 0 for validation dataset, 1 for raw_data_1
check_data=False        # check only first time you are using dataset


data_set=['validation_Q_l2','raw_data_1']

calculated_data_dir=path.join(path.dirname(__file__),'calculated data') 
data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),data_set[data_num])
calculated_data_file=path.join(calculated_data_dir,'vspace-'+data_set[data_num]+'.npy') 
data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)
u=data['u']
v=data['v']
w=data['w']


vspace1=np.zeros(np.shape(u))
vspace2=np.zeros(np.shape(u))
vspace3=np.zeros(np.shape(u))
delta=2.*math.pi/np.shape(u)[0]
x_max=np.shape(u)[0]-1
y_max=np.shape(u)[1]-1
z_max=np.shape(u)[2]-1
if n_elements>x_max:
    n_elements=x_max

def vel_der_ord2(vcomp,axis,p): # p is for point
    if axis=='x':
        s=np.array([1,0,0])  #s is for step
        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
    else: print ('wrong axis')   
    return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta

def vel_der_ord2new(vcomp,axis,p):
    if axis=='x':
        s=np.array([1,0,0])
        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[x_max,p[1]-s[1],p[2]-s[2]])/2./delta
        elif p[0]==x_max: return (vcomp[0,p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],y_max,p[2]-s[2]])/2./delta
        elif p[1]==y_max: return (vcomp[p[0]+s[0],0,p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],z_max])/2./delta
        elif p[2]==z_max: return (vcomp[p[0]+s[0],p[1]+s[1],0]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    else: print ('wrong axis')   
    return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta

def vel_der_ord4(vcomp,axis,p):
    if axis=='x':
        s=np.array([1,0,0])
        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
        elif p[0]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
        elif p[0]==x_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
        elif p[1]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
        elif p[1]==y_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
        elif p[2]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
        elif p[2]==z_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    else: print ('wrong axis')   
    return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
    
def vel_der_ord4new(vcomp,axis,p):
    if axis=='x':
        s=np.array([1,0,0])
        if p[0]==0: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[x_max,p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[x_max-1,p[1]-2*s[1],p[2]-2*s[2]])/12./delta
        elif p[0]==x_max: return (8*(vcomp[0,p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[1,p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
        elif p[0]==1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[x_max,p[1]-2*s[1],p[2]-2*s[2]])/12./delta
        elif p[0]==x_max-1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[0,p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],y_max,p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],y_max-1,p[2]-2*s[2]])/12./delta
        elif p[1]==y_max: return (8*(vcomp[p[0]+s[0],0,p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],1,p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
        elif p[1]==1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],y_max,p[2]-2*s[2]])/12./delta
        elif p[1]==y_max-1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],0,p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],z_max])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],z_max-1])/12./delta
        elif p[2]==z_max: return (8*(vcomp[p[0]+s[0],p[1]+s[1],0]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],1]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
        elif p[2]==1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],z_max])/12./delta
        elif p[2]==z_max-1: return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],0]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
    else: print ('wrong axis')   
    return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta

def vel_der_ord6new(vcomp,axis,p):
    if axis=='x':
        s=np.array([1,0,0])
        if p[0]==0: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[x_max,p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[x_max-1,p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[n_elements-3,p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[0]==x_max: return (45*(vcomp[0,p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[1,p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[2,p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[0]==1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[x_max,p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[x_max-1,p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[0]==x_max-1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[0,p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[1,p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[0]==2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[x_max,p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[0]==x_max-2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[0,p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],y_max,p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],y_max-1,p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],n_elements-3,p[2]-3*s[2]])/60./delta
        elif p[1]==y_max: return (45*(vcomp[p[0]+s[0],0,p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],1,p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],2,p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[1]==1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],y_max,p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],y_max-1,p[2]-3*s[2]])/60./delta
        elif p[1]==y_max-1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],0,p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],1,p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[1]==2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],y_max,p[2]-3*s[2]])/60./delta
        elif p[1]==y_max-2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],0,p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],z_max])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],z_max-1])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],n_elements-3])/60./delta
        elif p[2]==z_max: return (45*(vcomp[p[0]+s[0],p[1]+s[1],0]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],1]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],2]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[2]==1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],z_max])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],z_max-1])/60./delta
        elif p[2]==z_max-1: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],0]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],1]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
        elif p[2]==2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],z_max])/60./delta
        elif p[2]==z_max-2: return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],0]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta
    else: print ('wrong axis')   
    return (45*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],p[2]+3*s[2]]-vcomp[p[0]-3*s[0],p[1]-3*s[1],p[2]-3*s[2]])/60./delta


    

    
def D_matrix(point):
    if order_der_method==4:
        D=np.array([[vel_der_ord4(u,'x',point), vel_der_ord4(u,'y',point), vel_der_ord4(u,'z',point)],\
                    [vel_der_ord4(v,'x',point), vel_der_ord4(v,'y',point), vel_der_ord4(v,'z',point)],\
                    [vel_der_ord4(w,'x',point), vel_der_ord4(w,'y',point), vel_der_ord4(w,'z',point)]])
    elif order_der_method==2:
         D=np.array([[vel_der_ord2(u,'x',point), vel_der_ord2(u,'y',point), vel_der_ord2(u,'z',point)],\
                    [vel_der_ord2(v,'x',point), vel_der_ord2(v,'y',point), vel_der_ord2(v,'z',point)],\
                    [vel_der_ord2(w,'x',point), vel_der_ord2(w,'y',point), vel_der_ord2(w,'z',point)]])
    elif order_der_method==3:
         D=np.array([[vel_der_ord2new(u,'x',point), vel_der_ord2new(u,'y',point), vel_der_ord2new(u,'z',point)],\
                    [vel_der_ord2new(v,'x',point), vel_der_ord2new(v,'y',point), vel_der_ord2new(v,'z',point)],\
                    [vel_der_ord2new(w,'x',point), vel_der_ord2new(w,'y',point), vel_der_ord2new(w,'z',point)]])
    elif order_der_method==5:
         D=np.array([[vel_der_ord4new(u,'x',point), vel_der_ord4new(u,'y',point), vel_der_ord4new(u,'z',point)],\
                    [vel_der_ord4new(v,'x',point), vel_der_ord4new(v,'y',point), vel_der_ord4new(v,'z',point)],\
                    [vel_der_ord4new(w,'x',point), vel_der_ord4new(w,'y',point), vel_der_ord4new(w,'z',point)]])
    elif order_der_method==6:
         D=np.array([[vel_der_ord6new(u,'x',point), vel_der_ord6new(u,'y',point), vel_der_ord6new(u,'z',point)],\
                    [vel_der_ord6new(v,'x',point), vel_der_ord6new(v,'y',point), vel_der_ord6new(v,'z',point)],\
                    [vel_der_ord6new(w,'x',point), vel_der_ord6new(w,'y',point), vel_der_ord6new(w,'z',point)]])
    return D
    
def S_matrix(Dmatrix):
    return (Dmatrix+np.transpose(Dmatrix))/2.    
        
#   O is Omega matrix       
def O_matrix(Dmatrix):
    return (Dmatrix-np.transpose(Dmatrix))/2.    
    
def A_matrix(matS,matO):
    return np.dot(matS,matS)+np.dot(matO,matO)
    
def norm(matrix):
    mat=np.dot(matrix,np.transpose(matrix))
    return (mat[0,0]+mat[1,1]+mat[2,2])**0.5
    
def Q(normO,normS):
    return 0.5*(normO**2-normS**2)
    
def calc_Q(point):
    D=D_matrix(point)
    return Q(norm(O_matrix(D)),norm(S_matrix(D)))

def Lambda2(point):
    w, v = np.linalg.eigh(A_matrix(S_matrix(D_matrix(point)),O_matrix(D_matrix(point))))
    return w[1]    
    

print ('start calc')
stop1 = time.clock()
if to_calc_Q:
        
        
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=3
                matrix1=D_matrix(np.array([i,j,k]))
                vspace1[i,j,k]=calc_Q(np.array([i,j,k]))
    stop2=time.clock()            
    print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation, 2nd order')            
    
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=5
                matrix2=D_matrix(np.array([i,j,k]))
                vspace2[i,j,k]=calc_Q(np.array([i,j,k]))
    stop3=time.clock()            
    print ('\n',int((time.clock()-stop2)*10000)/10000.,'sec  Q criterion calculation,4th order')
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=6
                matrix3=D_matrix(np.array([i,j,k]))
                vspace3[i,j,k]=calc_Q(np.array([i,j,k]))
    stop4=time.clock()            
                 
    print ('\n',int((time.clock()-stop3)*10000)/10000.,'sec  Q criterion calculation, 6th order')
   
        
        
elif to_calc_Lambda2:
    
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=3
                matrix1=D_matrix(np.array([i,j,k]))
                vspace1[i,j,k]=Lambda2(np.array([i,j,k]))
    print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation, 2nd order')
                
    stop2=time.clock()  
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=5
                matrix2=D_matrix(np.array([i,j,k]))
                vspace2[i,j,k]=Lambda2(np.array([i,j,k]))
    print ('\n',int((time.clock()-stop2)*10000)/10000.,'sec  Q criterion calculation,4th order')
    stop3=time.clock()  
    for i in range(n_elements):
        for j in range(n_elements):
            for k in range(n_elements):
                order_der_method=6
                matrix3=D_matrix(np.array([i,j,k]))
                vspace3[i,j,k]=Lambda2(np.array([i,j,k]))
    stop4=time.clock() 
    print ('\n',int((time.clock()-stop3)*10000)/10000.,'sec  Q criterion calculation, 6th order')
    
    
    
    
    
calc_time=int((time.clock()-stop1)*10000)/10000.



abs24=abs(vspace2-vspace1)
abs26=abs(vspace3-vspace1)
abs46=abs(vspace3-vspace2)
print (np.sum(abs(matrix1)),np.sum(abs(matrix2)),np.sum(abs(matrix3)))
print (np.sum(abs24)/n_elements**3,np.sum(abs26)/n_elements**3, np.sum(abs46)/n_elements**3)   
print(np.sum(abs(vspace2))/n_elements**3,np.sum(abs(vspace3))/n_elements**3,np.sum(abs(vspace1))/n_elements**3)