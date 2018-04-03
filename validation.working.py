# -*- coding: utf-8 -*-
"""
Created on Sun Jun 04 18:31:55 2017

@author: Lisa
"""

import numpy as np
from numpy import *
import numpy as np
from math import *
import matplotlib.pyplot as plt

def vel_der_ord2(vcomp,axis,p): #p is for point
    if axis=='x':
        s=np.array([1,0,0])  #s is for s
        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1]] - vcomp[p[0],p[1]])/delta
        elif p[0]==x_max: return (vcomp[p[0],p[1]] - vcomp[p[0]-s[0],p[1]-s[1]])/delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
        elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
    else: print ('wrong axis')   
    return (vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])/2./delta

def vel_der_ord2new(vcomp,axis,p): #p is for point
    if axis=='x':
        s=np.array([1,0,0])  #s is for s
        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[x_max,p[1]-s[1]])/2./delta
        elif p[0]==x_max: return (vcomp[0,p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])/2./delta
    elif axis=='y':
        s=np.array([0,1,0])
        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],y_max,p[2]-s[2]])/2./delta
        elif p[1]==y_max: return (vcomp[p[0]+s[0],0,p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    elif axis=='z':
        s=np.array([0,0,1])
        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],z_max])/2./delta
        elif p[2]==z_max: return (vcomp[p[0]+s[0],p[1]+s[1],0]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
    else: print ('wrong axis')   
    return (vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])/2./delta
#
#def vel_der_ord4(vcomp,axis,p): #p is for point
#    if axis=='x':
#        s=np.array([1,0,0])  #s is for step
#        if p[0]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
#        elif p[0]==x_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
#        elif p[0]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#        elif p[0]==x_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#    elif axis=='y':
#        s=np.array([0,1,0])
#        if p[1]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
#        elif p[1]==y_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
#        elif p[1]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#        elif p[1]==y_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#    elif axis=='z':
#        s=np.array([0,0,1])
#        if p[2]==0: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]] - vcomp[p[0],p[1],p[2]])/delta
#        elif p[2]==z_max: return (vcomp[p[0],p[1],p[2]] - vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/delta
#        elif p[2]==1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#        elif p[2]==z_max-1: return (vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])/2./delta
#    else: print ('wrong axis')   
#    return (8*(vcomp[p[0]+s[0],p[1]+s[1],p[2]+s[2]]-vcomp[p[0]-s[0],p[1]-s[1],p[2]-s[2]])-vcomp[p[0]+2*s[0],p[1]+2*s[1],p[2]+2*s[2]]+vcomp[p[0]-2*s[0],p[1]-2*s[1],p[2]-2*s[2]])/12./delta
#    
def vel_der_ord4new(vcomp,axis,p): #p is for point
    if axis=='x':
        s=np.array([1,0,0])  #s is for step
        if p[0]==0: return (8*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[x_max,p[1]-s[1]])-vcomp[p[0]+2*s[0],p[1]+2*s[1]]+vcomp[(x_max-1),p[1]-2*s[1]])/12./delta
        elif p[0]==x_max: return (8*(vcomp[0,p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-vcomp[1,p[1]+2*s[1]]+vcomp[p[0]-2*s[0],p[1]-2*s[1]])/12./delta
        elif p[0]==1: return (8*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-vcomp[p[0]+2*s[0],p[1]+2*s[1]]+vcomp[x_max,p[1]-2*s[1]])/12./delta
        elif p[0]==x_max-1: return (8*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-vcomp[0,p[1]+2*s[1]]+vcomp[p[0]-2*s[0],p[1]-2*s[1]])/12./delta
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
    return (8*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-vcomp[p[0]+2*s[0],p[1]+2*s[1]]+vcomp[p[0]-2*s[0],p[1]-2*s[1]])/12./delta
#
#
def vel_der_ord6new(vcomp,axis,p): #p is for point
    if axis=='x':
        s=np.array([1,0,0])  #s is for step
        if p[0]==0: return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[x_max,p[1]-s[1]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1]]-vcomp[x_max-1,p[1]-2*s[1]])+vcomp[p[0]+3*s[0],p[1]+3*s[1]]-vcomp[x_max-2,p[1]-3*s[1]])/60./delta
        elif p[0]==x_max: return (45*(vcomp[0,p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[1,p[1]+2*s[1]]-vcomp[p[0]-2*s[0],p[1]-2*s[1]])+vcomp[2,p[1]+3*s[1]]-vcomp[p[0]-3*s[0],p[1]-3*s[1]])/60./delta
        elif p[0]==1: return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1]]-vcomp[x_max,p[1]-2*s[1],])+vcomp[p[0]+3*s[0],p[1]+3*s[1]]-vcomp[x_max-1,p[1]-3*s[1]])/60./delta
        elif p[0]==x_max-1: return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[0,p[1]+2*s[1]]-vcomp[p[0]-2*s[0],p[1]-2*s[1],])+vcomp[1,p[1]+3*s[1]]-vcomp[p[0]-3*s[0],p[1]-3*s[1]])/60./delta
        elif p[0]==2: return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1],]-vcomp[p[0]-2*s[0],p[1]-2*s[1]])+vcomp[p[0]+3*s[0],p[1]+3*s[1]]-vcomp[x_max,p[1]-3*s[1]])/60./delta
        elif p[0]==x_max-2: return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1]]-vcomp[p[0]-2*s[0],p[1]-2*s[1]])+vcomp[0,p[1]+3*s[1]]-vcomp[p[0]-3*s[0],p[1]-3*s[1]])/60./delta
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
    return (45*(vcomp[p[0]+s[0],p[1]+s[1]]-vcomp[p[0]-s[0],p[1]-s[1]])-9*(vcomp[p[0]+2*s[0],p[1]+2*s[1]]-vcomp[p[0]-2*s[0],p[1]-2*s[1]])+vcomp[p[0]+3*s[0],p[1]+3*s[1],]-vcomp[p[0]-3*s[0],p[1]-3*s[1]])/60./delta

M=500
errorlist1=[]
errorlist2=[]
errorlist3=[]
steplist=[]

for N in range(6,M+1):
    delta=2*pi/N
    u1=np.empty([N,3])
    u2=np.empty([N,3])
    
    for n in range(N):
        u1[n,0]=n*delta
        u2[n,0]=n*delta
        u2[n,1]=pi/2
        u2[n,2]=pi/2
        
    u_sin=np.sin(u1) 
    realder=np.cos(u2)    
    
    approximant2=np.empty([N,3])
    approximant4=np.empty([N,3])
    approximant6=np.empty([N,3])
    delta=2*pi/np.shape(u1)[0]
    steplist.append(delta)
    x_max=np.shape(u1)[0]-1

    
    for i in range(np.shape(u1)[0]):
        
        b=vel_der_ord2new(u_sin,'x',([i,0,0]))
        
        approximant2[i,0]=b
        
        c=vel_der_ord4new(u_sin,'x',([i,0,0]))
        approximant4[i,0]=c
        
        d=vel_der_ord6new(u_sin,'x',([i,0,0]))
        approximant6[i,0]=d
        
        
    error2ord=abs(realder-approximant2)
    error4ord=abs(realder-approximant4)
    error6ord=abs(realder-approximant6)
    errorlist1.append(error2ord[0,0])
    errorlist2.append(error4ord[0,0])
    errorlist3.append(error6ord[0,0])
    
y1=[]
y2=[]
y3=[]
#for x in steplist:  
#    y1.append(x**(-2))
#    y2.append(x**(-4))
#    y3.append(x**(-6))

for i in range(len(steplist)):
    steplist[i]=log(steplist[i])
    errorlist1[i]=log(errorlist1[i])
    errorlist2[i]=log(errorlist2[i])
    errorlist3[i]=log(errorlist3[i])
#    y1[i]=log(y1[i])
#    y2[i]=log(y2[i])
#    y3[i]=log(y3[i])
coeff1=np.polyfit(steplist,errorlist1,1)
coeff2=np.polyfit(steplist,errorlist2,1)
coeff3=np.polyfit(steplist,errorlist3,1)
print (coeff1,coeff2,coeff3)

plt.subplot(131)
plt.plot(steplist,errorlist1)
plt.ylim(-25,0)
plt.title('2nd order')

plt.subplot(132)
plt.plot(steplist,errorlist2)
plt.ylim(-25,0)
plt.title('4th order')

plt.subplot(133)
plt.plot(steplist,errorlist3)
plt.ylim(-25,0)
plt.title("6th order")
plt.show()

