#----------------------------------------------------------Modules used for plotting 
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.plotly as py
import plotly
from plotly.graph_objs import *
import plotly.figure_factory
from skimage import measure
##from pyevtk.hl import gridToVTK
#import matplotlib.pyplot
#graveyard pf unused module for now
#from mpl_toolkits.mplot3d import Axes3D


#plotly authentification, I can give you the access to the account just ask
plotly.tools.set_credentials_file(username='hunter139', api_key='lgN7Sd8dqPktT2wwpfCc')


q_threshold=0.16          # threshold for marching cubes algorithm 



xx=np.zeros((5,5,5))
un=xx.astype(dtype=str)
for i in range(5):
    for j in range(5):
        for k in range(5):
            un[i,j,k]=str(str(i)+str(j)+str(k))
print(un)




nx=2
ny=5
nz=8
mat=np.zeros((nx,ny,nz))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            mat[i,j,k]=1*i+2*j+4*k  
      
mat1x=np.concatenate((mat,np.zeros((2,np.shape(mat)[1],np.shape(mat)[2]))),axis=0)            
#print(mat1x)            
mat2x=np.concatenate((np.zeros((2,np.shape(mat)[1],np.shape(mat)[2])),mat),axis=0)           
derx=(mat1x-mat2x)/2.
derx[1,:,:]=mat1x[1,:,:]-mat1x[0,:,:]
derx[-2,:,:]=mat2x[-1,:,:]-mat2x[-2,:,:]
finx=derx[1:-1,:,:]
print(finx)

mat1y=np.concatenate((mat,np.zeros((np.shape(mat)[0],2,np.shape(mat)[2]))),axis=1)            
#print(mat1x)            
mat2y=np.concatenate((np.zeros((np.shape(mat)[0],2,np.shape(mat)[2])),mat),axis=1)           
dery=(mat1y-mat2y)/2.
dery[:,1,:]=mat1y[:,1,:]-mat1y[:,0,:]
dery[:,-2,:]=mat2y[:,-1,:]-mat2y[:,-2,:]
finy=dery[:,1:-1,:]
print(finy)

mat1z=np.concatenate((mat,np.zeros((np.shape(mat)[0],np.shape(mat)[1],2))),axis=2)                     
mat2z=np.concatenate((np.zeros((np.shape(mat)[0],np.shape(mat)[1],2)),mat),axis=2)          
derz=(mat1z-mat2z)/2.
derz[:,:,1]=mat1z[:,:,1]-mat1z[:,:,0]
derz[:,:,-2]=mat2z[:,:,-1]-mat2z[:,:,-2]
print(derz[:,:,1:-1])





def D_matrix2(point):
    return(np.array([[vel_der_ord2x(u,point), vel_der_ord2y(u,point), vel_der_ord2z(u,point)],\
                    [vel_der_ord2x(v,point), vel_der_ord2y(v,point), vel_der_ord2z(v,point)],\
                    [vel_der_ord2x(w,point), vel_der_ord2y(w,point), vel_der_ord2z(w,point)]]))
def D_matrix4(point):
    return(np.array([[vel_der_ord4x(u,point), vel_der_ord4y(u,point), vel_der_ord4z(u,point)],\
                    [vel_der_ord4x(v,point), vel_der_ord4y(v,point), vel_der_ord4z(v,point)],\
                    [vel_der_ord4x(w,point), vel_der_ord4y(w,point), vel_der_ord4z(w,point)]]))
    
    
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

def S_matrixold(Dmatrix):
    return (Dmatrix+np.transpose(Dmatrix))/2. 
def O_matrixold(Dmatrix):
    return (Dmatrix-np.transpose(Dmatrix))/2. 

def normnew(m):
    return (np.sum(m*m))**0.5
norm=normold #old is better :/
#normold(O_matrix(D_matrix([10,10,10]))) 

def Qnew(normO,normS):
    return 0.5*(normO*normO-normS*normS)
Q=Qold #old is better

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


def D_matrix2loop(point):
    return(np.array([[vel_der_ord2loopx(u,point), vel_der_ord2loopy(u,point), vel_der_ord2loopz(u,point)],\
                    [vel_der_ord2loopx(v,point), vel_der_ord2loopy(v,point), vel_der_ord2loopz(v,point)],\
                    [vel_der_ord2loopx(w,point), vel_der_ord2loopy(w,point), vel_der_ord2loopz(w,point)]]))
    
def D_matrix2(point):
    return(np.array([[vel_der_ord2x(u,point), vel_der_ord2y(u,point), vel_der_ord2z(u,point)],\
                    [vel_der_ord2x(v,point), vel_der_ord2y(v,point), vel_der_ord2z(v,point)],\
                    [vel_der_ord2x(w,point), vel_der_ord2y(w,point), vel_der_ord2z(w,point)]]))
    
def D_matrix4loop(point):
    return(np.array([[vel_der_ord4loopx(u,point), vel_der_ord4loopy(u,point), vel_der_ord4loopz(u,point)],\
                    [vel_der_ord4loopx(v,point), vel_der_ord4loopy(v,point), vel_der_ord4loopz(v,point)],\
                    [vel_der_ord4loopx(w,point), vel_der_ord4loopy(w,point), vel_der_ord4loopz(w,point)]]))
    
def D_matrix4(point):
    return(np.array([[vel_der_ord4x(u,point), vel_der_ord4y(u,point), vel_der_ord4z(u,point)],\
                    [vel_der_ord4x(v,point), vel_der_ord4y(v,point), vel_der_ord4z(v,point)],\
                    [vel_der_ord4x(w,point), vel_der_ord4y(w,point), vel_der_ord4z(w,point)]]))

def D_matrix6loop(point):
    return(np.array([[vel_der_ord6loopx(u,point), vel_der_ord6loopy(u,point), vel_der_ord6loopz(u,point)],\
                    [vel_der_ord6loopx(v,point), vel_der_ord6loopy(v,point), vel_der_ord6loopz(v,point)],\
                    [vel_der_ord6loopx(w,point), vel_der_ord6loopy(w,point), vel_der_ord6loopz(w,point)]]))

def D_matrix6(point):
    return(np.array([[vel_der_ord6x(u,point), vel_der_ord6y(u,point), vel_der_ord6z(u,point)],\
                    [vel_der_ord6x(v,point), vel_der_ord6y(v,point), vel_der_ord6z(v,point)],\
                    [vel_der_ord6x(w,point), vel_der_ord6y(w,point), vel_der_ord6z(w,point)]]))
    
#def vorticity(point):
#    if order_der_method==4:
#        i = vel_der_ord4y(w,point) - vel_der_ord4z(v,point)
#        j = -(vel_der_ord4x(w,point) - vel_der_ord4z(u,point))
#        k = vel_der_ord4x(v,point) - vel_der_ord4y(u,point)
#    elif order_der_method==2:
#        i = vel_der_ord2y(w,point) - vel_der_ord2z(v,point)
#        j = -(vel_der_ord2x(w,point) - vel_der_ord2z(u,point))
#        k = vel_der_ord2x(v,point) - vel_der_ord2y(u,point)
#    elif order_der_method==3:
#        i = vel_der_ord2loopy(w,point) - vel_der_ord2loopz(v,point)
#        j = -(vel_der_ord2loopx(w,point) - vel_der_ord2loopz(u,point))
#        k = vel_der_ord2loopx(v,point) - vel_der_ord2loopy(u,point)
#    elif order_der_method==5:
#        i = vel_der_ord4loopy(w,point) - vel_der_ord4loopz(v,point)
#        j = -(vel_der_ord4loopx(w,point) - vel_der_ord4loopz(u,point))
#        k = vel_der_ord4loopx(v,point) - vel_der_ord4loopy(u,point)
#    elif order_der_method==6:
#        i = vel_der_ord6loopy(w,point) - vel_der_ord6loopz(v,point)
#        j = -(vel_der_ord6loopx(w,point) - vel_der_ord6loopz(u,point))
#        k = vel_der_ord6loopx(v,point) - vel_der_ord6loopy(u,point)
#    strength = math.sqrt(i**2 + j**2 + k**2) 
#    return strength, i , j , k 
#
    
    if to_calc_Q:
        for i in range(3,n_elements+3):
            for j in range(3,n_elements+3):
                for k in range(3,n_elements+3):
                    vspace[i-3,j-3,k-3]=calc_Q(np.array([i,j,k]))
        print ('\n',int((time.clock()-stop1)*10000)/10000.,'sec  Q criterion calculation')
        highest_vorticity=np.amax(vspace)
    elif to_calc_Lambda2:
        for i in range(3,n_elements+3):
            for j in range(3,n_elements+3):
                for k in range(3,n_elements+3):
                    vspace[i-3,j-3,k-3]=Lambda2(np.array([i,j,k]))




















