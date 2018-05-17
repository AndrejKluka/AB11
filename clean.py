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

import math
stop=time.clock()
print ('\n',int((stop-start)*1000)/1000.,'sec -- imported modules')

''' To do list '''
#automate when we get more data

#---------------------------------------------------------General setup for program run
Visualization = True
to_save=True  
to_calc_Q=True       # if true will calc Q on cube with n_elements
to_calc_Lambda2=False   # if true will calc lambda2 on cube with n_elements
data_num=0            # 0 for validation dataset, 1 for raw_data_1, 2 for data_001,  3 for movie files
interval=110             # size of the cubes with which the program calculates Q/Lambda
frames=1              # frames to calc from movie
#65 -132sec
#100-125sec
#130-125sec
#160-114sec
#180-119sec
#256-154sec
optimal_intervals=[110,110,160,110]
interval=optimal_intervals[data_num]
data_set=['validation_Q_l2','raw_data_1','data_001','uvwp_00001.h5']

frame_names=['uvwp_00001.h5']
if data_num==3 and frames!=1:
    for i in range(frames):
        name='00'
        if (i+2)>=10:
            name='0'
        if (i+2)>=100:
            name=''
        name=name+str(i+2)
        frame_names.append('uvwp_00{}.h5' .format(name))
        
#   reading raw dataset and putting them into u,v,w arrays
data_set_file=path.join(path.join(path.dirname(__file__),'data sets'),data_set[data_num])  
movie_data=path.join(path.join(path.dirname(__file__),'data sets'),'Movie data')

#stuff just for fun
def print_statusline(msg: str):
    last_msg_length = len(print_statusline.last_msg) if hasattr(print_statusline, 'last_msg') else 0
    print(' ' * last_msg_length, end='\r')
    print(msg, end='\r')
    print_statusline.last_msg = msg
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
    s=np.zeros(np.shape(D)) #[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]+D[:,:,:,1,0])/2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]+D[:,:,:,2,0])/2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]+D[:,:,:,1,2])/2.
    s[:,:,:,0,0]=D[:,:,:,0,0]
    s[:,:,:,1,1]=D[:,:,:,1,1]
    s[:,:,:,2,2]=D[:,:,:,2,2]
    return(s)       

def O_matrix(D):
    s=np.zeros(np.shape(D)) #[0],np.shape(D)[1],np.shape(D)[2],3,3)
    s[:,:,:,0,1]=(D[:,:,:,0,1]-D[:,:,:,1,0])/2.
    s[:,:,:,1,0]=(D[:,:,:,0,1]-D[:,:,:,1,0])/-2.
    s[:,:,:,0,2]=(D[:,:,:,0,2]-D[:,:,:,2,0])/2.
    s[:,:,:,2,0]=(D[:,:,:,0,2]-D[:,:,:,2,0])/-2.
    s[:,:,:,2,1]=(D[:,:,:,2,1]-D[:,:,:,1,2])/2.
    s[:,:,:,1,2]=(D[:,:,:,2,1]-D[:,:,:,1,2])/-2.
    return(s) 

def norm_full(field):
    ft=np.transpose(field, (0,1,2,4,3))
    mat=np.matmul(field,ft)
    s=(mat[:,:,:,0,0]+mat[:,:,:,1,1]+mat[:,:,:,2,2])**0.5
    return s

def calc_Qfull(Dfield):
    qspace=0.5*(norm_full(O_matrix(Dfield))**2.-norm_full(S_matrix(Dfield))**2.)
    return qspace

def Lambda2full(Dfield):
    Ofield=O_matrix(Dfield)
    Sfield=S_matrix(Dfield)
    A=np.matmul(Sfield,Sfield)+np.matmul(Ofield,Ofield)
    lambda2_space=np.linalg.eigh(A)[0][:,:,:,1]
    return lambda2_space

times=1   
if data_num==3 and frames!=1:
    times=frames

for frame in range(times) :   
    if data_num==3:    
        filename =path.join(movie_data, frame_names[frame])
        f = h5py.File(filename, 'r')
        p=f[list(f.keys())[0]]
        u=f[list(f.keys())[1]]
        v=f[list(f.keys())[2]]
        w=f[list(f.keys())[3]]    
    else:    
        data=scipy.io.loadmat(data_set_file, mdict=None, appendmat=True)
        u=data['u']
        v=data['v']
        w=data['w']
    
    x_max=np.shape(u)[0]-1
    y_max=np.shape(u)[1]-1
    z_max=np.shape(u)[2]-1
    n_points=(x_max+1)*(y_max+1)*(z_max+1)
    points_calculated=0
    
    maxx=x_max
    if y_max>maxx:
        maxx=y_max
    elif z_max>maxx:
        maxx=z_max
    delta=2.*math.pi/(maxx+1)
    
    vspace=np.empty((u.shape),dtype='float32')
    vorticity_strength = np.empty((u.shape),dtype='float32')
    vorticity_x = np.empty((u.shape),dtype='float32')
    vorticity_y = np.empty((u.shape),dtype='float32')
    vorticity_z = np.empty((u.shape),dtype='float32')
    
    x=[0]
    y=[0]
    z=[0]
    axis_orig=[x,y,z]
    maxes=[x_max,y_max,z_max]
    interval=100
    for i in range(3):
        while axis_orig[i][-1]<maxes[i]:
            axis_orig[i].append(axis_orig[i][-1]+interval)
        axis_orig[i][-1]=maxes[i]
        if axis_orig[i][-1]-axis_orig[i][-2]<=20:
            axis_orig[i][-2]=axis_orig[i][-2]-30
    
    if to_calc_Q:
        method_of_choice=calc_Qfull
    elif to_calc_Lambda2:
        method_of_choice=Lambda2full
        
        
    print('start of calc')
    print_statusline(str(int(points_calculated/n_points*100))+'%')
    stop1 = time.clock()
    
    for i in range(len(axis_orig[0])-1):
        for j in range(len(axis_orig[1])-1):
            for k in range(len(axis_orig[2])-1):
                axis = copy.deepcopy(axis_orig)
                nah=[i,j,k]
                start=[0,0,0]
                end=[axis_orig[0][i+1]-axis_orig[0][i], axis_orig[1][j+1]-axis_orig[1][j], axis_orig[2][k+1]-axis_orig[2][k]]
                for xx in range(3):
                    if not axis[xx][nah[xx]]==0:
                        axis[xx][nah[xx]]-=3
                        start[xx]=3
                        end[xx]+=3
                    if not axis[xx][nah[xx]+1]==maxes[xx]:
                        axis[xx][nah[xx]+1]+=3
                
                Dfields=full_D_matrix( u[axis[0][i]:axis[0][i+1], axis[1][j]:axis[1][j+1], axis[2][k]:axis[2][k+1]], v[axis[0][i]:axis[0][i+1], axis[1][j]:axis[1][j+1], axis[2][k]:axis[2][k+1]],w[axis[0][i]:axis[0][i+1], axis[1][j]:axis[1][j+1], axis[2][k]:axis[2][k+1]],6)
                vorticity_strength[axis_orig[0][i]:axis_orig[0][i+1], axis_orig[1][j]:axis_orig[1][j+1], axis_orig[2][k]:axis_orig[2][k+1]] = Dfields[1][start[0]:end[0], start[1]:end[1], start[2]:end[2]]
                vorticity_x[axis_orig[0][i]:axis_orig[0][i+1], axis_orig[1][j]:axis_orig[1][j+1], axis_orig[2][k]:axis_orig[2][k+1]] = Dfields[2][start[0]:end[0], start[1]:end[1], start[2]:end[2]]
                vorticity_y[axis_orig[0][i]:axis_orig[0][i+1], axis_orig[1][j]:axis_orig[1][j+1], axis_orig[2][k]:axis_orig[2][k+1]] = Dfields[3][start[0]:end[0], start[1]:end[1], start[2]:end[2]]
                vorticity_z[axis_orig[0][i]:axis_orig[0][i+1], axis_orig[1][j]:axis_orig[1][j+1], axis_orig[2][k]:axis_orig[2][k+1]] = Dfields[4][start[0]:end[0], start[1]:end[1], start[2]:end[2]]
                zz=method_of_choice(Dfields[0])
                if to_calc_Q or to_calc_Lambda2:
                    vspace[axis_orig[0][i]:axis_orig[0][i+1], axis_orig[1][j]:axis_orig[1][j+1], axis_orig[2][k]:axis_orig[2][k+1]] = zz[start[0]:end[0], start[1]:end[1], start[2]:end[2]]
                points_calculated+=(-start[0]+end[0])*(-start[1]+end[1])*(-start[2]+end[2])
                print_statusline(str(int(points_calculated/n_points*100))+'%')        
     
    vorticity_x=np.nan_to_num(vorticity_x/vorticity_strength)
    vorticity_y=np.nan_to_num(vorticity_y/vorticity_strength)
    vorticity_z=np.nan_to_num(vorticity_z/vorticity_strength)
    print_statusline(str(100)+'%')
    print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  calculations done') 
    
    
    if to_save:
        if to_calc_Q: method='Q'
        elif to_calc_Lambda2: method='Lambda2'     
        xvtk = np.arange(0, vspace.shape[0])
        yvtk = np.arange(0, vspace.shape[1])
        zvtk = np.arange(0, vspace.shape[2])
        if data_num==3 and frames!=1:
            addon=frame_names[frame]
        else:
            addon=data_set[data_num]
        gridToVTK("./calculated data/" + addon + "-"+ method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_strength, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })
        #gridToVTK("C:\\Users\\Public\\Calculated_data\\" + data_set[data_num] + "-"+ method, xvtk, yvtk, zvtk, pointData = {method: vspace, "Vorticity normal": vorticity_strength, "Vorticity x" : vorticity_x , "Vorticity y" : vorticity_y , "Vorticity z" : vorticity_z })
        
    
    if Visualization : 
        os.chdir("C:\\Program Files\\ParaView 5.5.0-RC3-Qt5-Windows-64bit\\bin\\")
        os.system("pvpython.exe C:\\Users\\Public\\pv1.py")

