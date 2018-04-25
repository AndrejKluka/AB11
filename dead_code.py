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























