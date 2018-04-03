import time
start = time.clock()
#----------------------------------------------------------Modules used for plotting 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.plotly as py
import plotly
from plotly.graph_objs import *
import plotly.figure_factory
from pyevtk.hl import gridToVTK
#import matplotlib.pyplot
#graveyard pf unused module for now
#from mpl_toolkits.mplot3d import Axes3D


#plotly authentification, I can give you the access to the account just ask
plotly.tools.set_credentials_file(username='hunter139', api_key='lgN7Sd8dqPktT2wwpfCc')

#----------------------------------------------------------Modules for general life
from os import path
import scipy.io
import scipy
import numpy as np
from skimage import measure
import math
stop=time.clock()
print ('\n',int((stop-start)*1000)/1000.,'sec -- imported modules')




#---------------------------------------------------------General setup for program run
to_load=False          # if true will load already the last calculated Q or lambda dataset
to_plotly=False        # if true will send the plot to plotly website
to_matplot=False        # if true will use matplotlib to plot
n_elements=192        # number of elements on each side of cube calculated
to_calc_Q=False          # if true will calc Q on cube with n_elements
to_calc_Lambda2=True   # if true will calc lambda2 on cube with n_elements
to_calc_vorticity = True  #if true calculate vorticity
q_threshold=0.16          # threshold for marching cubes algorithm 
order_der_method=3     # only 2 or 4 are implemented 3 is 2 but new
data_num=1              # 0 for validation dataset, 1 for raw_data_1
check_data=False        # check only first time you are using dataset



data_set=['validation_Q_l2','raw_data_1']

#   reading raw dataset and putting them into u,v,w arrays
calculated_data_dir=path.join(path.dirname(__file__),'calculated data') 
data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),data_set[data_num])
calculated_data_file=path.join(calculated_data_dir,'vspace-'+data_set[data_num]+'.npy') 
data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)
u=data['u']
v=data['v']
w=data['w']

ushape = u.shape


if check_data:
    ok=True
    if not np.shape(u)==np.shape(v)==np.shape(w):
        ok=False
    if not np.shape(u)[0]==np.shape(u)[1]==np.shape(u)[2]:
        ok=False
    for i in range(np.shape(u)[0]):
        for j in range(np.shape(u)[0]):
            for k in range(np.shape(u)[0]):
                if math.isnan(u[i,j,k]) or math.isnan(v[i,j,k]) or math.isnan(w[i,j,k]):
                    print(i,j,k)
                    ok=False
    if ok:print('data is rectangular and ok')
    else:print('data not fine')

#vspace=np.zeros(np.shape(u))
vspace=np.zeros((n_elements,n_elements,n_elements))
vorticity_space = np.zeros((n_elements,n_elements,n_elements))
delta=2.*math.pi/np.shape(u)[0]
x_max=np.shape(u)[0]-1
y_max=np.shape(u)[1]-1
z_max=np.shape(u)[2]-1
if n_elements>x_max:
    n_elements=x_max

#   Definitions  for calculations
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
#   velocity gradient matrix
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

def vorticity(point):
    if order_der_method==4:
        i = vel_der_ord4(w,'y',point) - vel_der_ord4(v,'z',point)
        j = vel_der_ord4(w,'x',point) - vel_der_ord4(u,'z',point)
        k = vel_der_ord4(v,'x',point) - vel_der_ord4(u,'y',point)
    elif order_der_method==2:
        i = vel_der_ord2(w,'y',point) - vel_der_ord2(v,'z',point)
        j = vel_der_ord2(w,'x',point) - vel_der_ord2(u,'z',point)
        k = vel_der_ord2(v,'x',point) - vel_der_ord2(u,'y',point)
    elif order_der_method==3:
        i = vel_der_ord2new(w,'y',point) - vel_der_ord2new(v,'z',point)
        j = vel_der_ord2new(w,'x',point) - vel_der_ord2new(u,'z',point)
        k = vel_der_ord2new(v,'x',point) - vel_der_ord2new(u,'y',point)
    elif order_der_method==5:
        i = vel_der_ord4new(w,'y',point) - vel_der_ord4new(v,'z',point)
        j = vel_der_ord4new(w,'x',point) - vel_der_ord4new(u,'z',point)
        k = vel_der_ord4new(v,'x',point) - vel_der_ord4new(u,'y',point)
    elif order_der_method==6:
        i = vel_der_ord6new(w,'y',point) - vel_der_ord6new(v,'z',point)
        j = vel_der_ord6new(w,'x',point) - vel_der_ord6new(u,'z',point)
        k = vel_der_ord6new(v,'x',point) - vel_der_ord6new(u,'y',point)
    strength = math.sqrt(i**2 + j**2 + k**2) 
    return strength

if to_load:
    stop1 = time.clock()
    vspace=np.load(calculated_data_file)
    calc_time=int((time.clock()-stop1)*10000)/10000.
    print ('\n',calc_time,'sec  loaded calculation')
    highest_vorticity=np.amax(vspace) # need to be careful, here I assume that I load calculated Q values and no L2
else:
    print ('start calc')
    stop1 = time.clock()
    if to_calc_Q:
        for i in range(n_elements):
            for j in range(n_elements):
                for k in range(n_elements):
                    vspace[i,j,k]=calc_Q(np.array([i,j,k]))
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
        highest_vorticity=np.amax(vspace)
    elif to_calc_Lambda2:
        for i in range(n_elements):
            for j in range(n_elements):
                for k in range(n_elements):
                    vspace[i,j,k]=Lambda2(np.array([i,j,k]))
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Lambda2 calculation')
        highest_vorticity=np.amin(vspace)
    
    calc_time=int((time.clock()-stop1)*10000)/10000.
    np.save(calculated_data_file,vspace)  
if to_calc_vorticity:
    stop2 = time.clock()
    for i in range(n_elements):
            for j in range(n_elements):
                for k in range(n_elements):
                    vorticity_space[i,j,k]=vorticity(np.array([i,j,k]))
    print ('\n',int((time.clock()-stop2)*10000)/10000.,'sec  vorticity strength calculation')
    highest_vorticity=np.amax(vorticity_space)

 
   
    
#   plotting in plotly if set true      
if to_plotly:
    verts, simplices = measure.marching_cubes_classic(vspace, np.amax(vspace)*q_threshold)
    x,y,z = zip(*verts)
    #colormap=['rgb(255,105,180)','rgb(255,255,51)','rgb(0,191,255)']
    fig = plotly.figure_factory.create_trisurf(x=x, y=y, z=z, plot_edges=False,
                        #colormap=colormap,
                        simplices=simplices,
                        title="Vortex field")
    py.iplot(fig)
#   plotting in matplotlib if set true 
if to_matplot:
    stop1 = time.clock()        
    verts, faces = measure.marching_cubes_classic(vspace, np.amax(vspace)*q_threshold)
    stop2 = time.clock()   
    print ('\n',int((stop2-stop1)*10000)/10000.,'sec  marching cubes')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('g')
    ax.add_collection3d(mesh)
    #ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
                    #linewidth=0.1,antialiased=True)
                    #ax.set_edgecolor((255,255,255))
                    #cmap='Spectral'
    ax.set_xlim(0, n_elements)
    ax.set_ylim(0, n_elements)
    ax.set_zlim(0, n_elements)
    plt.show()




#   Saving calculation times to calctimes.txt
if not to_load:
    if to_calc_Q: method='Q       '
    elif to_calc_Lambda2: method='Lambda2 '
    wri=str('points ='+str(n_elements**3)+'   order of method='+str(order_der_method)+'  method='+method+'  time taken='+str(calc_time)+'sec    time per point='+str(calc_time/(n_elements**3.)*1000000)+'^10-6  ' )
    f=open('calctimes.txt','a')
    f.write(wri)
    f.write('\n')
    f.close()



vspace_shape = np.shape(vspace)      
xvtk = np.arange(0, vspace_shape[0])
yvtk = np.arange(0, vspace_shape[1])
zvtk = np.arange(0, vspace_shape[2])



gridToVTK("./calculated data/" + data_set[0] + "-" + str(n_elements) + "of" + str(np.shape(u)[0]) + "-" + method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity": vorticity_space})
