import time
start = time.clock()


#----------------------------------------------------------Modules for general life
from os import path
import scipy.io
import scipy
import numpy as np
#from pyevtk.hl import gridToVTK

import math
stop=time.clock()
print ('\n',int((stop-start)*1000)/1000.,'sec -- imported modules')


#---------------------------------------------------------General setup for program run
to_load=False          # if true will load already the last calculated Q or lambda dataset
to_save=False
<<<<<<< HEAD
to_plotly=False        # if true will send the plot to plotly website
to_matplot=False        # if true will use matplotlib to plot

<<<<<<< HEAD
n_elements_x=95       #in x-direction
n_elements_y=95        # in y-direction
n_elements_z=95   # number of elements on each side of cube calculated- z-direction
=======


n_elements_x=255       #in x-direction
n_elements_y=255        # in y-direction
n_elements_z=767   # number of elements on each side of cube calculated- z-direction
>>>>>>> 2d731150056e12b6ef7164f0055e573834f74e14
=======
n_elements_x=100       #in x-direction
n_elements_y=100        # in y-direction
n_elements_z=100   # number of elements on each side of cube calculated- z-direction
>>>>>>> edb4c20bb3e675176679807d9b0084fe13a56193
q_threshold=0.16          # threshold for marching cubes algorithm     
<<<<<<< HEAD
to_calc_Q=True        # if true will calc Q on cube with n_elements
=======
<<<<<<< HEAD
to_calc_Q=True        # if true will calc Q on cube with n_elements
=======
to_calc_Q=False        # if true will calc Q on cube with n_elements
>>>>>>> 5a19d77e6eb8da8f93b7de45cc0438ba0fe31dcb
>>>>>>> e05b1f21863a6e3bc3a297c9819be3906c036702
to_calc_Lambda2=False   # if true will calc lambda2 on cube with n_elements
to_calc_vorticity = True  #if true calculate vorticity
order_der_method=6     #2,4,6 are with looping in 2,4,6 orders respectetively
to_loop=False         # True if the data loops 
<<<<<<< HEAD
<<<<<<< HEAD
data_num=0              # 0 for validation dataset, 1 for raw_data_1, 2 for data_001
=======
data_num=2              # 0 for validation dataset, 1 for raw_data_1, 2 for data_001
>>>>>>> 2d731150056e12b6ef7164f0055e573834f74e14
=======
to_save_calctime=False  #saving calc times to text file
data_num=1              # 0 for validation dataset, 1 for raw_data_1, 2 for data_001
>>>>>>> edb4c20bb3e675176679807d9b0084fe13a56193

check_data=False        # check only first time you are using dataset
 

data_set=['validation_Q_l2','raw_data_1','data_001']

#   reading raw dataset and putting them into u,v,w arrays
calculated_data_dir=path.join(path.dirname(__file__),'calculated data') 
data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),data_set[data_num])
calculated_data_file=path.join(calculated_data_dir,'vspace-'+data_set[data_num]+'.npy')
'''calculated_vorticity_file=path.join(calculated_data_dir,'vorticity-'+data_set[data_num]+'.npy') 
calculated_vorticityx_file=path.join(calculated_data_dir,'vorticity-x-'+data_set[data_num]+'.npy') 
calculated_vorticityy_file=path.join(calculated_data_dir,'vorticity-y-'+data_set[data_num]+'.npy')
calculated_vorticityz_file=path.join(calculated_data_dir,'vorticity-z-'+data_set[data_num]+'.npy') ''' 
 
data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)
u=data['u']
v=data['v']
w=data['w']



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
    if ok:print('data is cubical and ok')
    else:print('data not fine')

vspace=np.zeros(np.shape(u))
x_max=np.shape(u)[0]-1
y_max=np.shape(u)[1]-1
z_max=np.shape(u)[2]-1


maxx=x_max
if y_max>maxx:
    maxx=y_max
elif z_max>maxx:
    maxx=z_max
delta=2.*math.pi/maxx

   
if n_elements_x>x_max:
    n_elements_x=x_max
elif n_elements_y>y_max:
    n_elements_y=y_max
elif n_elements_z>z_max:
    n_elements_z=z_max

vspace=np.zeros((n_elements_x,n_elements_y,n_elements_z))
vorticity_space = np.zeros((n_elements_x,n_elements_y,n_elements_z))
vorticity_x = np.zeros((n_elements_x,n_elements_y,n_elements_z))
vorticity_y = np.zeros((n_elements_x,n_elements_y,n_elements_z))
vorticity_z = np.zeros((n_elements_x,n_elements_y,n_elements_z))

#stuff for finite diff method
def ord2_full_mat(mat):
    derx=np.zeros((np.shape(mat)))
    derx[1:-1,:,:]=(mat[2:,:,:]-mat[:-2,:,:])/2
    derx[0,:,:]=mat[1,:,:]-mat[0,:,:]
    derx[-1,:,:]=mat[-1,:,:]-mat[-2,:,:]
    
    dery=np.zeros((np.shape(mat)))
    dery[:,1:-1,:]=(mat[:,2:,:]-mat[:,:-2,:])/2
    dery[:,0,:]=mat[:,1,:]-mat[:,0,:]
    dery[:,-1,:]=mat[:,-1,:]-mat[:,-2,:]
    
    derz=np.zeros((np.shape(mat)))
    derz[:,:,1:-1]=(mat[:,:,2:]-mat[:,:,:-2])/2
    derz[:,:,0]=mat[:,:,1]-mat[:,:,0]
    derz[:,:,-1]=mat[:,:,-1]-mat[:,:,-2]
    return derx/delta, dery/delta, derz/delta

def ord4_full_mat(mat):
    derx=np.zeros((np.shape(mat)))
    derx[2:-2,:,:]=(8*(mat[3:-1,:,:]-mat[1:-3,:,:])-mat[4:,:,:]+mat[:-4,:,:])/12
    derx[0,:,:]=mat[1,:,:]-mat[0,:,:]
    derx[1,:,:]=(mat[2,:,:]-mat[0,:,:])/2
    derx[-2,:,:]=(mat[-1,:,:]-mat[-3,:,:])/2
    derx[-1,:,:]=mat[-1,:,:]-mat[-2,:,:]
    
    dery=np.zeros((np.shape(mat)))
    dery[:,2:-2,:]=(8*(mat[:,3:-1,:]-mat[:,1:-3,:])-mat[:,4:,:]+mat[:,:-4,:])/12
    dery[:,0,:]=mat[:,1,:]-mat[:,0,:]
    dery[:,1,:]=(mat[:,2,:]-mat[:,0,:])/2
    dery[:,-2,:]=(mat[:,-1,:]-mat[:,-3,:])/2
    dery[:,-1,:]=mat[:,-1,:]-mat[:,-2,:]
    
    derz=np.zeros((np.shape(mat)))
    derz[:,:,2:-2]=(8*(mat[:,:,3:-1]-mat[:,:,1:-3])-mat[:,:,4:]+mat[:,:,:-4])/12
    derz[:,:,0]=mat[:,:,1]-mat[:,:,0]
    derz[:,:,1]=(mat[:,:,2]-mat[:,:,0])/2
    derz[:,:,-2]=(mat[:,:,-1]-mat[:,:,-3])/2
    derz[:,:,-1]=mat[:,:,-1]-mat[:,:,-2]
    
    return derx/delta, dery/delta, derz/delta

def ord6_full_mat(mat):
    derx=np.zeros((np.shape(mat)))
    derx[3:-3,:,:]=(45*(mat[4:-2,:,:]-mat[2:-4,:,:])-9*(mat[5:-1,:,:]-mat[1:-5,:,:])+mat[6:,:,:]-mat[:-6,:,:])/60
    derx[0,:,:]=mat[1,:,:]-mat[0,:,:]
    derx[1,:,:]=(mat[2,:,:]-mat[0,:,:])/2
    derx[2,:,:]=(8*(mat[3,:,:]-mat[1,:,:])-mat[4,:,:]+mat[0,:,:])/12
    derx[-3,:,:]=(8*(mat[-2,:,:]-mat[-4,:,:])-mat[-1,:,:]+mat[-5,:,:])/12
    derx[-2,:,:]=(mat[-1,:,:]-mat[-3,:,:])/2
    derx[-1,:,:]=mat[-1,:,:]-mat[-2,:,:]
    
    dery=np.zeros((np.shape(mat)))
    dery[:,3:-3,:]=(45*(mat[:,4:-2,:]-mat[:,2:-4,:])-9*(mat[:,5:-1,:]-mat[:,1:-5,:])+mat[:,6:,:]-mat[:,:-6,:])/60
    dery[:,0,:]=mat[:,1,:]-mat[:,0,:]
    dery[:,1,:]=(mat[:,2,:]-mat[:,0,:])/2
    dery[:,2,:]=(8*(mat[:,3,:]-mat[:,1,:])-mat[:,4,:]+mat[:,0,:])/12
    dery[:,-3,:]=(8*(mat[:,-2,:]-mat[:,-4,:])-mat[:,-1,:]+mat[:,-5,:])/12
    dery[:,-2,:]=(mat[:,-1,:]-mat[:,-3,:])/2
    dery[:,-1,:]=mat[:,-1,:]-mat[:,-2,:]
    
    derz=np.zeros((np.shape(mat)))
    derz[:,:,3:-3]=(45*(mat[:,:,4:-2]-mat[:,:,2:-4])-9*(mat[:,:,5:-1]-mat[:,:,1:-5])+mat[:,:,6:]-mat[:,:,:-6])/60
    derz[:,:,0]=mat[:,:,1]-mat[:,:,0]
    derz[:,:,1]=(mat[:,:,2]-mat[:,:,0])/2
    derz[:,:,2]=(8*(mat[:,:,3]-mat[:,:,1])-mat[:,:,4]+mat[:,:,0])/12
    derz[:,:,-3]=(8*(mat[:,:,-2]-mat[:,:,-4])-mat[:,:,-1]+mat[:,:,-5])/12
    derz[:,:,-2]=(mat[:,:,-1]-mat[:,:,-3])/2
    derz[:,:,-1]=mat[:,:,-1]-mat[:,:,-2]
    
    return derx/delta, dery/delta, derz/delta



#stuff for concatenating matrixes
def full_D_matrix(u,v,w,order_of_method):
    if order_of_method==2:
        method=ord2_full_mat
    if order_of_method==4:
        method=ord4_full_mat
    if order_of_method==6:
        method=ord6_full_mat
    deru=method(u)
    derv=method(v)
    derw=method(w)
    f_line=np.concatenate((np.reshape(deru[0],(np.shape(deru[0])[0],np.shape(deru[0])[1],np.shape(deru[0])[2],1)), np.reshape(derv[0],(np.shape(derv[0])[0],np.shape(derv[0])[1],np.shape(derv[0])[2],1)), np.reshape(derw[0],(np.shape(derw[0])[0],np.shape(derw[0])[1],np.shape(derw[0])[2],1))),axis=3)
    s_line=np.concatenate((np.reshape(deru[1],(np.shape(deru[1])[0],np.shape(deru[1])[1],np.shape(deru[1])[2],1)), np.reshape(derv[1],(np.shape(derv[1])[0],np.shape(derv[1])[1],np.shape(derv[1])[2],1)), np.reshape(derw[1],(np.shape(derw[1])[0],np.shape(derw[1])[1],np.shape(derw[1])[2],1))),axis=3)
    t_line=np.concatenate((np.reshape(deru[2],(np.shape(deru[2])[0],np.shape(deru[2])[1],np.shape(deru[2])[2],1)), np.reshape(derv[2],(np.shape(derv[2])[0],np.shape(derv[2])[1],np.shape(derv[2])[2],1)), np.reshape(derw[2],(np.shape(derw[2])[0],np.shape(derw[2])[1],np.shape(derw[2])[2],1))),axis=3)
    i = derw[1] - derv[2]
    j = deru[2] - derw[0]
    k = derv[0] - deru[1]
    strength=(i**2. + j**2. + k**2.)**0.5
    gradient_tensor=np.concatenate(( np.reshape(f_line,(np.shape(f_line)[0],np.shape(f_line)[1],np.shape(f_line)[2],3,1)),  np.reshape(s_line,(np.shape(s_line)[0],np.shape(s_line)[1],np.shape(s_line)[2],3,1)),  np.reshape(t_line,(np.shape(t_line)[0],np.shape(t_line)[1],np.shape(t_line)[2],3,1))), axis=4)
    return gradient_tensor, strength, i ,j ,k


#extending velocity fields in all directions if the data repeats 
def extend_matrix(matrix):
    zzx=np.concatenate((np.array(matrix[(np.shape(matrix)[0]-3):np.shape(matrix)[0],:,:]),matrix[:,:,:],np.array(matrix[0:3,:,:])), axis=0)
    matrix=zzx
    zzx=np.concatenate((np.array(matrix[:,(np.shape(matrix)[1]-3):np.shape(matrix)[1],:]),matrix[:,:,:],np.array(matrix[:,0:3,:])), axis=1)
    matrix=zzx
    zzx=np.concatenate((np.array(matrix[:,:,(np.shape(matrix)[2]-3):np.shape(matrix)[2]]),matrix[:,:,:],np.array(matrix[:,:,0:3])), axis=2)
    return(zzx)
 
if to_loop:
    u=extend_matrix(u)
    v=extend_matrix(v)
    w=extend_matrix(w)
    
    
#   Definitions  for calculations
def vel_der_ord2loopx(vcomp,p):
    return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
def vel_der_ord2loopy(vcomp,p):
    return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
def vel_der_ord2loopz(vcomp,p):
    return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta

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


def vel_der_ord4loopx(vcomp,p):
    return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
def vel_der_ord4loopy(vcomp,p):
    return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
def vel_der_ord4loopz(vcomp,p):
    return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],p[2]-2])/12./delta

def vel_der_ord4x(vcomp,p):
    if p[0]==0: return (vcomp[p[0]+1,p[1],p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-1,p[1],p[2]])/delta
    elif p[0]==1 or p[0]==x_max-1: return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
    return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
def vel_der_ord4y(vcomp,p):
    if p[1]==0: return (vcomp[p[0],p[1]+1,p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1]-1,p[2]])/delta
    elif p[1]==1 or p[1]==y_max-1: return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
    return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
def vel_der_ord4z(vcomp,p):
    if p[2]==0: return (vcomp[p[0],p[1],p[2]+1] - vcomp[p[0],p[1],p[2]])/delta
    elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1],p[2]-1])/delta
    elif p[2]==1 or p[2]==z_max-1: return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta
    return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],p[2]-2])/12./delta
    

def vel_der_ord6loopx(vcomp,p):
    return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta
def vel_der_ord6loopy(vcomp,p):
    return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta
def vel_der_ord6loopz(vcomp,p):
    return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],p[2]-3])/60./delta

def vel_der_ord6x(vcomp,p):
    if p[0]==0: return (vcomp[p[0]+1,p[1],p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-1,p[1],p[2]])/delta
    elif p[0]==1 or p[0]==x_max-1: return (vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])/2./delta
    elif p[0]==2 or p[0]==x_max-2: return (8*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-vcomp[p[0]+2,p[1],p[2]]+vcomp[p[0]-2,p[1],p[2]])/12./delta
    return (45*(vcomp[p[0]+1,p[1],p[2]]-vcomp[p[0]-1,p[1],p[2]])-9*(vcomp[p[0]+2,p[1],p[2]]-vcomp[p[0]-2,p[1],p[2]])+vcomp[p[0]+3,p[1],p[2]]-vcomp[p[0]-3,p[1],p[2]])/60./delta
def vel_der_ord6y(vcomp,p):
    if p[1]==0: return (vcomp[p[0],p[1]+1,p[2]] - vcomp[p[0],p[1],p[2]])/delta
    elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1]-1,p[2]])/delta
    elif p[1]==1 or p[1]==y_max-1: return (vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])/2./delta
    elif p[1]==2 or p[1]==y_max-2: return (8*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-vcomp[p[0],p[1]+2,p[2]]+vcomp[p[0],p[1]-2,p[2]])/12./delta
    return (45*(vcomp[p[0],p[1]+1,p[2]]-vcomp[p[0],p[1]-1,p[2]])-9*(vcomp[p[0],p[1]+2,p[2]]-vcomp[p[0],p[1]-2,p[2]])+vcomp[p[0],p[1]+3,p[2]]-vcomp[p[0],p[1]-3,p[2]])/60./delta
def vel_der_ord6z(vcomp,p):
    if p[2]==0: return (vcomp[p[0],p[1],p[2]+1] - vcomp[p[0],p[1],p[2]])/delta
    elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0],p[1],p[2]-1])/delta
    elif p[2]==1 or p[2]==z_max-1: return (vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])/2./delta
    elif p[2]==2 or p[2]==z_max-2: return (8*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-vcomp[p[0],p[1],p[2]+2]+vcomp[p[0],p[1],p[2]-2])/12./delta
    return (45*(vcomp[p[0],p[1],p[2]+1]-vcomp[p[0],p[1],p[2]-1])-9*(vcomp[p[0],p[1],p[2]+2]-vcomp[p[0],p[1],p[2]-2])+vcomp[p[0],p[1],p[2]+3]-vcomp[p[0],p[1],p[2]-3])/60./delta

#   velocity gradient matrix   
def D_matrix2(point):
    xu=vel_der_ord2x(u,point)
    xv=vel_der_ord2x(v,point)
    xw=vel_der_ord2x(w,point)
    yu=vel_der_ord2y(u,point)
    yv=vel_der_ord2y(v,point)
    yw=vel_der_ord2y(w,point)
    zu=vel_der_ord2z(u,point)    
    zv=vel_der_ord2z(v,point)
    zw=vel_der_ord2z(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k
    
def D_matrix2loop(point):
    xu=vel_der_ord2loopx(u,point)
    xv=vel_der_ord2loopx(v,point)
    xw=vel_der_ord2loopx(w,point)
    yu=vel_der_ord2loopy(u,point)
    yv=vel_der_ord2loopy(v,point)
    yw=vel_der_ord2loopy(w,point)
    zu=vel_der_ord2loopz(u,point)    
    zv=vel_der_ord2loopz(v,point)
    zw=vel_der_ord2loopz(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k
   
def D_matrix4(point):
    xu=vel_der_ord4x(u,point)
    xv=vel_der_ord4x(v,point)
    xw=vel_der_ord4x(w,point)
    yu=vel_der_ord4y(u,point)
    yv=vel_der_ord4y(v,point)
    yw=vel_der_ord4y(w,point)
    zu=vel_der_ord4z(u,point)    
    zv=vel_der_ord4z(v,point)
    zw=vel_der_ord4z(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k
    
def D_matrix4loop(point):
    xu=vel_der_ord4loopx(u,point)
    xv=vel_der_ord4loopx(v,point)
    xw=vel_der_ord4loopx(w,point)
    yu=vel_der_ord4loopy(u,point)
    yv=vel_der_ord4loopy(v,point)
    yw=vel_der_ord4loopy(w,point)
    zu=vel_der_ord4loopz(u,point)    
    zv=vel_der_ord4loopz(v,point)
    zw=vel_der_ord4loopz(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k
    

def D_matrix6loop(point):
    xu=vel_der_ord6loopx(u,point)
    xv=vel_der_ord6loopx(v,point)
    xw=vel_der_ord6loopx(w,point)
    yu=vel_der_ord6loopy(u,point)
    yv=vel_der_ord6loopy(v,point)
    yw=vel_der_ord6loopy(w,point)
    zu=vel_der_ord6loopz(u,point)    
    zv=vel_der_ord6loopz(v,point)
    zw=vel_der_ord6loopz(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k
   
def D_matrix6(point):    
    xu=vel_der_ord6x(u,point)
    xv=vel_der_ord6x(v,point)
    xw=vel_der_ord6x(w,point)
    yu=vel_der_ord6y(u,point)
    yv=vel_der_ord6y(v,point)
    yw=vel_der_ord6y(w,point)
    zu=vel_der_ord6z(u,point)    
    zv=vel_der_ord6z(v,point)
    zw=vel_der_ord6z(w,point)
    i=yw-zv
    j=-xw+zu
    k=xv-yu
    strength= math.sqrt(i**2 + j**2 + k**2)
    return(np.array([[xu, yu,zu],\
                    [xv, yv, zv],\
                    [xw,yw , zw]])), strength, i, j, k


if to_loop:        
    if order_der_method==2:   D_matrix=D_matrix2loop 
    elif order_der_method==4:    D_matrix=D_matrix4loop   
    elif order_der_method==6:    D_matrix=D_matrix6loop
else: 
    if order_der_method==2:   D_matrix=D_matrix2
    elif order_der_method==4:    D_matrix=D_matrix4   
    elif order_der_method==6:    D_matrix=D_matrix6


def S_matrix(D):
    s=np.zeros(np.shape(D))#[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    #print (D)
    return(s)       

def O_matrix(D):
    s=np.zeros(np.shape(D))#[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]-D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]-D[:,:,:,1,0])/2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]-D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]-D[:,:,:,2,0])/2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]-D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]-D[:,:,:,1,2])/2.
    #print(s)
    return(s)        

   
def S_matrixold(D):  
    s=np.zeros((3,3))
    s[0,1]=s[1,0]=(D[0,1]+D[1,0])/2.
    s[0,2]=s[2,0]=(D[0,2]+D[2,0])/2.
    s[2,1]=s[1,2]=(D[2,1]+D[1,2])/2.
    return(s)       
 
#   O is Omega matrix        
def O_matrixold(D):
    s=np.zeros((3,3))
    s[0,1]=s[1,0]=(D[0,1]-D[1,0])/2.
    s[0,2]=s[2,0]=(D[0,2]-D[2,0])/2.
    s[2,1]=s[1,2]=(D[2,1]-D[1,2])/2.
    return(s)  
 
def A_matrix(matS,matO):
    return np.dot(matS,matS)+np.dot(matO,matO)

def norm(m):
    mat=np.dot(m,np.transpose(m))
    return (mat[0,0]+mat[1,1]+mat[2,2])**0.5

def norm_full(field):
    ft=np.transpose(field, (0,1,2,4,3))
    mat=np.matmul(field,ft)
    s=(mat[:,:,:,0,0]+mat[:,:,:,1,1]+mat[:,:,:,2,2])**0.5
    return s
   
def Q(normO,normS):
    return 0.5*(normO**2.-normS**2.)
  
def calc_Q(point):
    D=D_matrix(point) 
    return Q(norm(O_matrixold(D[0])),norm(S_matrixold(D[0]))),D[1], D[2], D[3], D[4]   #q value, vorticity strenght, vorticity i,j,k

def calc_Qfull(Dfield):
    qspace=0.5*(norm_full(O_matrix(Dfield[0]))**2.-norm_full(S_matrix(Dfield[0]))**2.)
    return qspace, Dfield[1], Dfield[2], Dfield[3], Dfield[4]
    
 
    
def calc_Qx(D):
    return Q(norm(O_matrixold(D)),norm(S_matrixold(D)))  #q value, vorticity strenght, vorticity i,j,k


def Lambda2(point):
    D=D_matrix(point)
    w, v = np.linalg.eigh(A_matrix(S_matrixold(D[0]),O_matrixold(D[0])))
    return w[1], D[1], D[2], D[3], D[4]

def Lambda2full(Dfield):
    Ofield=O_matrix(Dfield[0])
    Sfield=S_matrix(Dfield[0])
    A=np.matmul(Sfield,Sfield)+np.matmul(Ofield,Ofield)
    lambda2_space=np.linalg.eigh(A)[0][:,:,:,1]
    return lambda2_space, Dfield[1], Dfield[2], Dfield[3], Dfield[4]


def Lambda2x(D):
    w, v = np.linalg.eigh(A_matrix(S_matrixold(D),O_matrixold(D)))
    return w[1]



if to_load:
    stop1 = time.clock()
    vspace=np.load(calculated_data_file)
    '''vorticity_space=np.load(calculated_vorticity_file)
    vorticity_x=np.load(calculated_vorticityx_file)
    vorticity_y=np.load(calculated_vorticityy_file)
    vorticity_z=np.load(calculated_vorticityz_file)'''
    calc_time=int((time.clock()-stop1)*10000)/10000.
   # print ('\n',calc_time,'sec  loaded calculation')
    highest_vorticity=np.amax(vspace) # need to be careful, here I assume that I load calculated Q values and no L2
else:
    #print ('start calc')
    stop1 = time.clock()
    if to_calc_Q and to_calc_vorticity:
        for i in range(0,x_max):
            for j in range(0,y_max):
                for k in range(0,z_max):
                    Qandvorticity=calc_Q(np.array([i,j,k]))
                    vspace[i,j,k]=Qandvorticity[0]
                    vorticity_space[i,j,k]=Qandvorticity[1]
                    vorticity_x[i,j,k]=Qandvorticity[2]
                    vorticity_y[i,j,k]=Qandvorticity[3]
                    vorticity_z[i,j,k]=Qandvorticity[4]
<<<<<<< HEAD
                    if (i,j,k) == (10,10,10):
                        print (vspace, Qandvorticity[0])
                        break
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
        
=======
                    if i==10 and j==10 and k==10: print (vspace[10,10,10])
<<<<<<< HEAD
        #print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
>>>>>>> 2d731150056e12b6ef7164f0055e573834f74e14
=======
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
>>>>>>> edb4c20bb3e675176679807d9b0084fe13a56193
        highest_vorticity=np.amax(vspace)
    elif to_calc_Q:
        for i in range(n_elements_x):           
            for j in range(n_elements_y):
                for k in range(n_elements_z):
                    Qandvorticity=calc_Q(np.array([i,j,k]))
                    vspace[i,j,k]=Qandvorticity[0]
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
        highest_vorticity=np.amax(vspace)        
    elif to_calc_Lambda2 and to_calc_vorticity: 
        for i in range(n_elements_x):
            for j in range(n_elements_y):
                for k in range(n_elements_z):
                    Lambdaandvorticity=Lambda2(np.array([i,j,k]))
                    vspace[i,j,k]=Lambdaandvorticity[0]
                    vorticity_space[i,j,k]=Lambdaandvorticity[1]
                    vorticity_x[i,j,k]=Lambdaandvorticity[2]
                    vorticity_y[i,j,k]=Lambdaandvorticity[3]
                    vorticity_z[i,j,k]=Lambdaandvorticity[4]
                    
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Lambda2 calculation')
        highest_vorticity=np.amin(vspace)
    elif to_calc_Lambda2:
        for i in range(n_elements_x):
             for j in range(n_elements_y):
                 for k in range(n_elements_z):
                     Lambdaandvorticity=Lambda2(np.array([i,j,k]))
                     vspace[i,j,k]=Lambdaandvorticity[0]
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Lambda2 calculation')
        highest_vorticity=np.amin(vspace)
    calc_time=int((time.clock()-stop1)*10000)/10000.
    

if to_save: np.save(calculated_data_file,vspace)  



'''
n_elements=30
vspace=np.zeros((n_elements,n_elements,n_elements))
print('okee')

stop1 = time.clock()    
jaa=full_D_matrix(u[0:n_elements,0:n_elements,0:n_elements],v[0:n_elements,0:n_elements,0:n_elements],w[0:n_elements,0:n_elements,0:n_elements],6)
zz=calc_Qfull(jaa)
print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  new D')
print(np.shape(zz))
#print (jaa)

#yy=Lambda2new(jaa)

stop1 = time.clock()
for i in range(n_elements):
    for j in range(n_elements):
        for k in range(n_elements):  
            vspace[i,j,k]=calc_Qx(D_matrix6([i,j,k]))
print(np.shape(vspace))
print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  old D')
print(np.sum(vspace[:-3,:-3,:-3]-zz[:-3,:-3,:-3]))
'''






#   Saving calculation times to calctimes.txt
<<<<<<< HEAD
<<<<<<< HEAD
#if not to_load:
#    if to_calc_Q: method='Q'
#    elif to_calc_Lambda2: method='Lambda2'
#    wri=str('points ='+str(n_elements**3)+'   order of method='+str(order_der_method)+'  method='+method+'  time taken='+str(calc_time)+'sec    time per point='+str(calc_time/(n_elements**3.)*1000000)+'^10-6  ' )
#    f=open('calctimes.txt','a')
#    f.write(wri)
#    f.write('\n')
#    f.close()


#vspace_shape = np.shape(vspace)      
#xvtk = np.arange(0, vspace_shape[0])
#yvtk = np.arange(0, vspace_shape[1])
#zvtk = np.arange(0, vspace_shape[2])
#
=======
if not to_load:
=======
if to_save_calctime:
>>>>>>> edb4c20bb3e675176679807d9b0084fe13a56193
    if to_calc_Q: method='Q'
    elif to_calc_Lambda2: method='Lambda2'
    wri=str('points ='+str(n_elements_x*n_elements_y*n_elements_z)+'   order of method='+str(order_der_method)+'  method='+method+'  time taken='+str(calc_time)+'sec    time per point='+str(calc_time/(n_elements_x*n_elements_y*n_elements_z)*1000000)+'^10-6  ' )
    f=open('calctimes.txt','a')
    f.write(wri)
    f.write('\n')
    f.close()


vspace_shape = np.shape(vspace)      
xvtk = np.arange(0, vspace_shape[0])
yvtk = np.arange(0, vspace_shape[1])
zvtk = np.arange(0, vspace_shape[2])

<<<<<<< HEAD
>>>>>>> 2d731150056e12b6ef7164f0055e573834f74e14
#gridToVTK("./calculated data/" + data_set[data_num] + "-" + str(n_elements) + "of" + str(np.shape(u)[0]) + "-" + method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_space, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })
=======
gridToVTK("./calculated data/" + data_set[data_num] + "-" + str(n_elements) + "of" + str(np.shape(u)[0]) + "-" + method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_space, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })
>>>>>>> edb4c20bb3e675176679807d9b0084fe13a56193

