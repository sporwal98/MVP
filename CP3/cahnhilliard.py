import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.animation as anim

import pandas as pd
def main(phi0,dx,dt, size = 50,a = 0.1,K= 0.1,M=0.1, nsteps = 100000):
    r = size
    c = size
    phi = (np.random.random((size,size))*0.2-0.1)+phi0
    print(phi)
    mu = calc_mu(phi,dx,a,K) 

    #'''
    #Figure Setups
    fig = plt.figure()
    im = plt.imshow(phi, animated = True,vmin = -1,vmax = +1)
    plt.colorbar(im)
    #'''
    freeE = []
    t = []

    st = 0
    while(st<nsteps):
        #print(st)

        mu = calc_mu(phi,dx,a,K)
        phi+= dt*M*nabla2(mu,dx)
        if(np.mod(st,1000) == 0):
            #print(st)
            #'''
            #Animation
            f = open('mu.dat','w')
            for i in range(r):
                for j in range(c):
                    f.write('%d %d %lf\n'%(i,j,phi[i,j]))
            f.close()
            plt.cla()
            im = plt.imshow(phi,animated = True,vmin = -1, vmax = +1)
            plt.draw()
            plt.pause(0.001)
            #'''
            #CALCULATING FREE ENERGY
            f = get_f(phi,dx,a,K)
            freeE.append(f)
            t.append(st)
            
            #'''
        #print(phi)
            
        
        st+=1
        if(np.mod(st,10000)==0):
            print(st)
    '''
    plt.figure()
    plt.plot(t,freeE)
    plt.xlabel('Time')
    plt.ylabel('Free Energy')
    #plt.show()
    plt.savefig('FEvsT.png')
    #'''

    dic = {'Time':t,'FreeEnergy':freeE}
    df = pd.DataFrame.from_dict(dic)
    df.to_csv('che_freeEvsT_'+str(phi0)+'_dt'+str(dt)+'.csv')
    return

def nabla2(fn,dx):
    return (np.roll(fn,+1,axis = 0) + np.roll(fn,-1,axis = 0) + np.roll(fn,+1,axis = 1) + np.roll(fn,-1,axis = 1) - 4*fn)/dx**2

def calc_mu(phi,dx,a,K):
    return (-a*phi+a*phi**3-K*nabla2(phi,dx))

def get_f(phi,dx,a,K):
    return np.sum((-a*(phi**2)/2 + a*(phi**4)/4 +K*( der_x(phi,dx)**2+der_y(phi,dx)**2 )/2))


def der_x(fn,dx):
    return (np.roll(fn,1,axis = 0)-np.roll(fn,-1,axis =0))/(2*dx)

def der_y(fn,dx):
    return (np.roll(fn,1,axis = 1)-np.roll(fn,-1,axis =1))/(2*dx)
