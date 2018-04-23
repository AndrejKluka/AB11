import time
start = time.clock()
#----------------------------------------------------------Modules used for plotting 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.plotly as py
import plotly
from plotly.graph_objs import *
import plotly.figure_factory
##from pyevtk.hl import gridToVTK
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


# improve method 6 it does not work with 0,0,0

#---------------------------------------------------------General setup for program run
idk=False
to_load=False          # if true will load already the last calculated Q or lambda dataset
to_plotly=False        # if true will send the plot to plotly website
to_matplot=False        # if true will use matplotlib to plot
n_elements=10        # number of elements on each side of cube calculated
to_calc_Q=True          # if true will calc Q on cube with n_elements
to_calc_Lambda2=False   # if true will calc lambda2 on cube with n_elements
to_calc_vorticity = False  #if true calculate vorticity
q_threshold=0.16          # threshold for marching cubes algorithm 
order_der_method=5       # 2,4 are without looping, 3,5,6 are with looping in 2,4,6 orders respectetively
data_num=1              # 0 for validation dataset, 1 for raw_data_1
check_data=False        # check only first time you are using dataset

#if order_der_method==2 or order_der_method==4:loop=False    
loop=False
#elif order_der_method==3 or order_der_method==5 or order_der_method==6:loop=True

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
vorticity_x = np.zeros((n_elements,n_elements,n_elements))
vorticity_y = np.zeros((n_elements,n_elements,n_elements))
vorticity_z = np.zeros((n_elements,n_elements,n_elements))
delta=2.*math.pi/np.shape(u)[0]
x_max=np.shape(u)[0]-1
y_max=np.shape(u)[1]-1
z_max=np.shape(u)[2]-1
if n_elements>x_max:
    n_elements=x_max
#extending velocity fields in all directions if the data repeats 
def extend_matrix(matrix):
    zzx=np.concatenate((np.array(matrix[(np.shape(matrix)[0]-3):np.shape(matrix)[0],:,:]),matrix[:,:,:],np.array(matrix[0:3,:,:])), axis=0)
    matrix=zzx
    zzx=np.concatenate((np.array(matrix[:,(np.shape(matrix)[1]-3):np.shape(matrix)[1],:]),matrix[:,:,:],np.array(matrix[:,0:3,:])), axis=1)
    matrix=zzx
    zzx=np.concatenate((np.array(matrix[:,:,(np.shape(matrix)[2]-3):np.shape(matrix)[2]]),matrix[:,:,:],np.array(matrix[:,:,0:3])), axis=2)
    return(zzx)

if loop:
    u=extend_matrix(u)
    v=extend_matrix(v)
    w=extend_matrix(w)
    
    
    


#   Definitions  for calculations
def vel_der_ord2x(vcomp,p):
    if p[0]==0: return (vcomp[p[0]+1,p[1],p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-1,p[1],p[2]])/delta
    return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
def vel_der_ord2y(vcomp,p):
    if p[1]==0: return (vcomp[p[0],p[1]+1,p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1]-1,p[2]])/delta
    return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
def vel_der_ord2z(vcomp,p):
    if p[2]==0: return (vcomp[p[0],p[1],p[2]+1] - vcomp[p[0],p[1],p[2]])/delta
    elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1],p[2]-1])/delta  
    return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta

def vel_der_ord2loopx(vcomp,p):
    if p[0]==0: return (vcomp[p[0]+1,p[1],p[2]] - vcomp[x_max,p[1],p[2]])/delta
    elif p[0]==x_max: return (vcomp[0,p[1],p[2]] - vcomp[p[0]-1,p[1],p[2]])/delta
    return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
def vel_der_ord2loopy(vcomp,p):
    if p[1]==0: return (vcomp[p[0],p[1]+1,p[2]] - vcomp[y_max,p[1],p[2]])/delta
    elif p[1]==y_max: return (vcomp[p[0],0,p[2]] - vcomp[p[0],p[1]-1,p[2]])/delta
    return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
def vel_der_ord2loopz(vcomp,p):
    if p[2]==0: return (vcomp[p[0],p[1],p[2]+1] - vcomp[z_max,p[1],p[2]])/delta
    elif p[2]==z_max: return (vcomp[p[0],p[1],0] - vcomp[p[0],p[1],p[2]-1])/delta  
    return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta



def vel_der_ord4x(vcomp,p):
    if p[0]==0: return (vcomp[p[0]+1,p[1],p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-1,p[1],p[2]])/delta
    elif p[0]==1 or p[0]==x_max-1: return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
    return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
def vel_der_ord4y(vcomp,p):
    if p[1]==0: return (vcomp[p[0],p[1]+1,p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1]-1,p[2]])/delta
    elif p[1]==1 or p[0]==y_max-1: return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
    return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
def vel_der_ord4z(vcomp,p):
    if p[2]==0: return (vcomp[p[0],p[1],p[2]+1] - vcomp[p[0],p[1],p[2]])/delta
    elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1],p[2]-1])/delta
    elif p[2]==1 or p[0]==z_max-1: return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta
    return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],p[2]-2])/12./delta

def vel_der_ord4loopx(vcomp,p):
    if p[0]==0: return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[x_max,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[x_max-1,p[1],p[2]])/12./delta
    elif p[0]==x_max: return (8*(vcomp[0,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[1,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
    elif p[0]==1: return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[x_max,p[1],p[2]])/12./delta
    elif p[0]==x_max-1: return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[0,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
    return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
def vel_der_ord4loopy(vcomp,p):
    if p[1]==0: return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],y_max,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],y_max-1,p[2]])/12./delta
    elif p[1]==y_max: return (8*(vcomp[p[0],0,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],1,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
    elif p[1]==1: return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],y_max,p[2]])/12./delta
    elif p[1]==y_max-1: return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],0,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
    return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
def vel_der_ord4loopz(vcomp,p):
    if p[2]==0: return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],z_max])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],z_max-1])/12./delta
    elif p[2]==z_max: return (8*(vcomp[p[0],p[1],0]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],1]+vcomp[p[0],p[1],p[2]-2])/12./delta
    elif p[2]==1: return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],z_max])/12./delta
    elif p[2]==z_max-1: return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],0]+vcomp[p[0],p[1],p[2]-2])/12./delta
    return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],p[2]-2])/12./delta
    


def vel_der_ord6loopx(vcomp,p):
    if p[0]==0: return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[x_max,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[x_max-1,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[x_max-2,p[1],p[2]])/60./delta
    elif p[0]==x_max: return (45*(vcomp[0,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[1,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[2,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta
    elif p[0]==1: return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[x_max,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[x_max-1,p[1],p[2]])/60./delta
    elif p[0]==x_max-1: return(45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[0,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[1,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta
    elif p[0]==2: return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[x_max,p[1],p[2]])/60./delta
    elif p[0]==x_max-2: return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[0,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta
    return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta

def vel_der_ord6loopy(vcomp,p):
    if p[1]==0: return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],y_max,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],y_max-1,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],y_max-2,p[2]])/60./delta
    elif p[1]==y_max: return (45*(vcomp[p[0],0,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],1,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],2,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta
    elif p[1]==1: return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],y_max,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],y_max-1,p[2]])/60./delta
    elif p[1]==y_max-1: return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],0,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],1,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta
    elif p[1]==2: return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],y_max,p[2]])/60./delta
    elif p[1]==y_max-2: return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],0,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta
    return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta

def vel_der_ord6loopz(vcomp,p):
    if p[2]==0: return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],z_max])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],z_max-1])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],z_max-2])/60./delta
    elif p[2]==z_max: return (45*(vcomp[p[0],p[1],0]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],1]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],2]-vcomp[p[0],p[1],p[2]-3])/60./delta
    elif p[2]==1: return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],z_max])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],z_max-2])/60./delta
    elif p[2]==z_max-1: return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],1]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],2]-vcomp[p[0],p[1],p[2]-3])/60./delta
    elif p[2]==2: return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],z_max])/60./delta
    elif p[2]==z_max-2: return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],0]-vcomp[p[0],p[1],p[2]-3])/60./delta
    return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],p[2]-3])/60./delta


#   velocity gradient matrix
def D_matrix2(point):
    return(np.array([[vel_der_ord2x(u,point), vel_der_ord2y(u,point), vel_der_ord2z(u,point)],\
                    [vel_der_ord2x(v,point), vel_der_ord2y(v,point), vel_der_ord2z(v,point)],\
                    [vel_der_ord2x(w,point), vel_der_ord2y(w,point), vel_der_ord2z(w,point)]]))
def D_matrix2loop(point):
    return(np.array([[vel_der_ord2loopx(u,point), vel_der_ord2loopy(u,point), vel_der_ord2loopz(u,point)],\
                    [vel_der_ord2loopx(v,point), vel_der_ord2loopy(v,point), vel_der_ord2loopz(v,point)],\
                    [vel_der_ord2loopx(w,point), vel_der_ord2loopy(w,point), vel_der_ord2loopz(w,point)]]))
def D_matrix4(point):
    return(np.array([[vel_der_ord4x(u,point), vel_der_ord4y(u,point), vel_der_ord4z(u,point)],\
                    [vel_der_ord4x(v,point), vel_der_ord4y(v,point), vel_der_ord4z(v,point)],\
                    [vel_der_ord4x(w,point), vel_der_ord4y(w,point), vel_der_ord4z(w,point)]]))    
def D_matrix4loop(point):
    return(np.array([[vel_der_ord4loopx(u,point), vel_der_ord4loopy(u,point), vel_der_ord4loopz(u,point)],\
                    [vel_der_ord4loopx(v,point), vel_der_ord4loopy(v,point), vel_der_ord4loopz(v,point)],\
                    [vel_der_ord4loopx(w,point), vel_der_ord4loopy(w,point), vel_der_ord4loopz(w,point)]]))

def D_matrix6loop(point):
    return(np.array([[vel_der_ord6loopx(u,point), vel_der_ord6loopy(u,point), vel_der_ord6loopz(u,point)],\
                    [vel_der_ord6loopx(v,point), vel_der_ord6loopy(v,point), vel_der_ord6loopz(v,point)],\
                    [vel_der_ord6loopx(w,point), vel_der_ord6loopy(w,point), vel_der_ord6loopz(w,point)]]))

    
def vorticity(point):
    if order_der_method==5:
        i = vel_der_ord4y(w,point) - vel_der_ord4z(v,point)
        j = -(vel_der_ord4x(w,point) - vel_der_ord4z(u,point))
        k = vel_der_ord4x(v,point) - vel_der_ord4y(u,point)
#    elif order_der_method==2:
#        i = vel_der_ord2y(w,point) - vel_der_ord2z(v,point)
#        j = -(vel_der_ord2x(w,point) - vel_der_ord2z(u,point))
#        k = vel_der_ord2x(v,point) - vel_der_ord2y(u,point)
#    elif order_der_method==3:
#        i = vel_der_ord2loopy(w,point) - vel_der_ord2loopz(v,point)
#        j = -(vel_der_ord2loopx(w,point) - vel_der_ord2loopz(u,point))
#        k = vel_der_ord2loopx(v,point) - vel_der_ord2loopy(u,point)
#    elif order_der_method==4:
#        i = vel_der_ord4loopy(w,point) - vel_der_ord4loopz(v,point)
#        j = -(vel_der_ord4loopx(w,point) - vel_der_ord4loopz(u,point))
#        k = vel_der_ord4loopx(v,point) - vel_der_ord4loopy(u,point)
#    elif order_der_method==6:
#        i = vel_der_ord6loopy(w,point) - vel_der_ord6loopz(v,point)
#        j = -(vel_der_ord6loopx(w,point) - vel_der_ord6loopz(u,point))
#        k = vel_der_ord6loopx(v,point) - vel_der_ord6loopy(u,point)
    strength = math.sqrt(i**2 + j**2 + k**2) 
    return strength, i , j , k 


    
if order_der_method==2:   D_matrix=D_matrix2 
elif order_der_method==3:    D_matrix=D_matrix2loop
elif order_der_method==4:    D_matrix=D_matrix4
elif order_der_method==5:    D_matrix=D_matrix4loop   
elif order_der_method==6:    D_matrix=D_matrix6loop
 
   
def S_matrixold(Dmatrix):
    return (Dmatrix+np.transpose(Dmatrix))/2. 
def S_matrixnew(D):    
    D[0,1]=D[1,0]=(D[0,1]+D[1,0])/2.
    D[0,2]=D[2,0]=(D[0,2]+D[2,0])/2.
    D[2,1]=D[1,2]=(D[2,1]+D[1,2])/2.
    return(D)       
S_matrix=S_matrixnew #new is faster
 
   
#   O is Omega matrix       
def O_matrixold(Dmatrix):
    return (Dmatrix-np.transpose(Dmatrix))/2.  
def O_matrixnew(D):
    s=np.zeros((3,3))
    s[0,1]=s[1,0]=(D[0,1]-D[1,0])/2.
    s[0,2]=s[2,0]=(D[0,2]-D[2,0])/2.
    s[2,1]=s[1,2]=(D[2,1]-D[1,2])/2.
    return(s)       
O_matrix=O_matrixnew #new is faster  
 
   
def A_matrix(matS,matO):
    return np.dot(matS,matS)+np.dot(matO,matO)
    

def normold(m):
    mat=np.dot(m,np.transpose(m))
    return (mat[0,0]+mat[1,1]+mat[2,2])**0.5
def normnew(m):
    return (np.sum(m*m))**0.5
norm=normold #old is better :/
#normold(O_matrix(D_matrix([10,10,10]))) 

   
def Qold(normO,normS):
    return 0.5*(normO**2-normS**2)
def Qnew(normO,normS):
    return 0.5*(normO*normO-normS*normS)
Q=Qold #old is better
    
def calc_Q(point):
    D=D_matrix(point)
    return Q(norm(O_matrix(D)),norm(S_matrix(D)))

def Lambda2(point):
    w, v = np.linalg.eigh(A_matrix(S_matrix(D_matrix(point)),O_matrix(D_matrix(point))))
    return w[1]

def timed_calc_Q(point):
    stop0=time.clock()
    print ('\n',int((time.clock()-stop0)*10000000)/10.,'mikrosec  time check')
    stop1 = time.clock()
    D=D_matrix(point)
    print ('\n',int((time.clock()-stop1)*10000000)/10.,'mikrosec  D calculation')
    stop2 = time.clock()
    O=O_matrix(D)
    print ('\n',int((time.clock()-stop2)*10000000)/10.,'mikrosec  O calculation')
    stop3 = time.clock()
    S=S_matrix(D)
    print ('\n',int((time.clock()-stop3)*10000000)/10.,'mikrosec  S calculation')
    stop4 = time.clock() 
    N1=norm(O)
    print ('\n',int((time.clock()-stop4)*10000000)/10.,'mikrosec  Norm(O) calculation')
    stop5 = time.clock()
    N2=norm(S)
    print ('\n',int((time.clock()-stop5)*10000000)/10.,'mikrosec  Norm(S) calculation')
    stop6 = time.clock()
    QQ=Q(N1,N2)
    print ('\n',int((time.clock()-stop6)*10000000)/10.,'mikrosec  Q calculation')
    return QQ


#timed_calc_Q([10,10,10])
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
                    vorticity_space[i,j,k]=vorticity(np.array([i,j,k]))[0]
                    vorticity_x[i,j,k]=vorticity(np.array([i,j,k]))[1]
                    vorticity_y[i,j,k]=vorticity(np.array([i,j,k]))[2]
                    vorticity_z[i,j,k]=vorticity(np.array([i,j,k]))[3]
                    
                    
                    
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
    if to_calc_Q: method='Q'
    elif to_calc_Lambda2: method='Lambda2'
    wri=str('points ='+str(n_elements**3)+'   order of method='+str(order_der_method)+'  method='+method+'  time taken='+str(calc_time)+'sec    time per point='+str(calc_time/(n_elements**3.)*1000000)+'^10-6  ' )
    f=open('calctimes.txt','a')
    f.write(wri)
    f.write('\n')
    f.close()


if idk:
    vspace_shape = np.shape(vspace)      
    xvtk = np.arange(0, vspace_shape[0])
    yvtk = np.arange(0, vspace_shape[1])
    zvtk = np.arange(0, vspace_shape[2])


vspace_shape = np.shape(vspace)      
xvtk = np.arange(0, vspace_shape[0])
yvtk = np.arange(0, vspace_shape[1])
zvtk = np.arange(0, vspace_shape[2])

#gridToVTK("./calculated data/" + data_set[data_num] + "-" + str(n_elements) + "of" + str(np.shape(u)[0]) + "-" + method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_space, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })

