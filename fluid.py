# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 16:53:03 2020

@author: Janek
"""

import imageio
import numpy as np
from PIL import Image, ImageDraw 
import matplotlib
Nx=600
Ny=200
c=1
R=20
cx=150
cy=100
u0=0.05
epsilon=0.001
Ly=Ny-1
visco=u0*R/1000
tau=3*visco + 0.5
uin=np.array([u0,0])
e=np.array([[ 0, 0],[1,0],[ 0, 1],[-1,0],[0, -1],[1, 1],[-1, 1],[-1,-1],[1, -1]])
space=np.zeros([Nx,Ny])
f=np.ones([9,Nx,Ny])
w=np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
feq=np.copy(f)
fcol=np.copy(f)
u=np.zeros([2,Nx,Ny])
rho=np.ones([Nx,Ny])
y=np.arange(200)
filenames=[]
def macro_density():
    global f, rho
    rho=np.zeros([Nx,Ny])
    for a in range(9):
        rho+=f[a]
    #print("macro_density ", rho)
def macro_velocity():
    global u, rho
    u=np.zeros([2,Nx,Ny])
    for a in range(9):
        u[0]+=1/rho *f[a]*e[a][0]
        u[1]+=1/rho *f[a]*e[a][1]
def equilibrium_f(u):
    global feq
    #print("local: ", u[0][0][0], u[1][0][0] )
    for a in range(9):
        feq[a]=w[a]*rho*(1 + 3*(e[a][0]*u[0]+e[a][1]*u[1])/c**2 +9/2*((e[a][0]*u[0]+e[a][1]*u[1])/c**2)**2 - 3/2*(u[1]*u[1]+u[0]*u[0])/c**2)
def collision():
    global fcol, f
    print(tau)
    fcol = f - (f - feq)/tau
def streaming():
    global f
    for i in range(9):
        f[i,:,:] = np.roll(np.roll(fcol[i,:,:],e[i,0],axis=0),e[i,1],axis=1)

def boundary():
    global f,rho,u
   
    f[6][Nx-1,:]=f[6][Nx-2,:]
    f[3][Nx-1,:]=f[3][Nx-2,:]   #overflow on the right wall
    f[7][Nx-1,:]=f[7][Nx-2,:]
    macro_density()
    macro_velocity()
    u[0][0,:]=source()

    rho[0,:]=(2*(f[3][0]+f[6][0]+f[7][0])+f[0][0]+f[2][0]+f[4][0])/(1-np.sqrt(u[0][0,:]**2+u[1][0,:]**2))
    equilibrium_f(u)
    f[1][0,:]=feq[1][0,:]
    f[5][0,:]=feq[5][0,:]      #left wall 
    f[8][0,:]=feq[8][0,:]
    
    
    
    #f[6][:,Ny-1]=f[6][:,1]
    #f[2][:,Ny-1]=f[2][:,1]  #periodic top-bottom
    #f[5][:,Ny-1]=f[5][:,1]
    #f[7][:,1]=f[7][:,Ny-1]
    #f[4][:,1]=f[4][:,Ny-1]  #periodic top-bottom
    #f[8][:,1]=f[8][:,Ny-1]
    
   
    #print("boundary ",rho)
def cylinder_collision():
    global fcol,f
    reverse = [e.tolist().index((-e[i]).tolist()) for i in range(9)] 
    cylinder = np.fromfunction(lambda x,y: (x-cx)**2+(y-cy)**2<R**2, (Nx,Ny))
    for i in range(9):
        fcol[i,cylinder] = f[reverse[i],cylinder]
def display(t):
    matplotlib.pyplot.imshow(np.sqrt(u[0]**2+u[1]**2),cmap="plasma",norm=matplotlib.colors.Normalize(vmin=0, vmax=0.08, clip=False))
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.savefig("test{}.png".format(t))
    matplotlib.pyplot.close()
    filenames.append('test{}.png'.format(t))
def source():
    global u
    return(u0*(1+epsilon*np.sin(2*np.pi*y/Ly)))
   


u[0]=u0
equilibrium_f(u)
f=feq
#print(f[:,0,0])
for l in range(600):
    for t in range(100):
        boundary()
        collision()
        print(u[:,100,100])
        cylinder_collision()
        
        streaming()
    display(l)
    print(l)
 
images=[]
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('C:/Users/Janek/Desktop/ComputerModeling/Fluid_flow/3Re220.gif', images)

    

    