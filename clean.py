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

''' To do list '''
#see why lambda 2 is all positive
#Think about adding some smart loading
#add script for paraview thingy and automate
#make it calculate whole datasets


#---------------------------------------------------------General setup for program run
to_save=False
n_elements_x=100       #in x-direction
n_elements_y=100        # in y-direction
n_elements_z=100   # number of elements on each side of cube calculated- z-direction    
to_calc_Q=False        # if true will calc Q on cube with n_elements
to_calc_Lambda2=False   # if true will calc lambda2 on cube with n_elements
data_num=1              # 0 for validation dataset, 1 for raw_data_1, 2 for data_001


data_set=['validation_Q_l2','raw_data_1','data_001']

#   reading raw dataset and putting them into u,v,w arrays
data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),data_set[data_num]) 
data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)
u=data['u']
v=data['v']
w=data['w']

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
    
#calculating gradients with whole matrixes
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

def S_matrix(D):
    s=np.zeros(np.shape(D))#[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    return(s)       

def O_matrix(D):
    s=np.zeros(np.shape(D))#[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]-D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]-D[:,:,:,1,0])/2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]-D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]-D[:,:,:,2,0])/2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]-D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]-D[:,:,:,1,2])/2.
    return(s) 

def norm_full(field):
    ft=np.transpose(field, (0,1,2,4,3))
    mat=np.matmul(field,ft)
    s=(mat[:,:,:,0,0]+mat[:,:,:,1,1]+mat[:,:,:,2,2])**0.5
    return s

def calc_Qfull(Dfield):
    qspace=0.5*(norm_full(O_matrix(Dfield[0]))**2.-norm_full(S_matrix(Dfield[0]))**2.)
    return qspace, Dfield[1], Dfield[2], Dfield[3], Dfield[4]

def Lambda2full(Dfield):
    Ofield=O_matrix(Dfield[0])
    Sfield=S_matrix(Dfield[0])
    A=np.matmul(Sfield,Sfield)+np.matmul(Ofield,Ofield)
    lambda2_space=np.linalg.eigh(A)[0][:,:,:,1]
    return lambda2_space, Dfield[1], Dfield[2], Dfield[3], Dfield[4]
    


n_elements=130
vspace=np.zeros((n_elements,n_elements,n_elements))
print('okee')

stop1 = time.clock()    
jaa=full_D_matrix(u[0:n_elements,0:n_elements,0:n_elements],v[0:n_elements,0:n_elements,0:n_elements],w[0:n_elements,0:n_elements,0:n_elements],6)
zz=calc_Qfull(jaa)
print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  new D')
print(np.shape(zz))






vspace_shape = np.shape(vspace)      
xvtk = np.arange(0, vspace_shape[0])
yvtk = np.arange(0, vspace_shape[1])
zvtk = np.arange(0, vspace_shape[2])

#gridToVTK("./calculated data/" + data_set[data_num] + "-" + str(n_elements) + "of" + str(np.shape(u)[0]) + "-" + method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_space, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })






