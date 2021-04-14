import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

import matplotlib.animation as anim


def main(size = 50, nsteps = 5000, ini = 'R'):
    #Game of Life Simulation
    state = np.zeros((size,size))
    state = init_state(ini,size)

    (r,c) = state.shape

    xt = []
    yt = []
    comt = []

    #'''
    #Figure Setup
    fig = plt.figure()
    im = plt.imshow(state, animated = True)
    #'''
    eq_time = -1
    stepsfromeq = 0
    active = count_active(state)
    
    st = 0
    while(st<nsteps):
        #Time-step
        nn = aliveneighbours(state)

        xs,ys = np.where(state == 1)
        if((sum(xs == 0) == 0)  & (sum(xs == size-1)==0) & (sum(ys ==0)==0) & (sum(ys == size -1)==0)):
            #Checking if touching border of board
            com = com_g(xs,ys)
            xt.append(com[0])
            yt.append(com[1])
            comt.append(st)

        #Checking transitions
        lonely = np.logical_and(state == 1, nn<2)
        crowded = np.logical_and(state == 1, nn>3)
        alive = np.logical_and(state == 1, np.logical_or(nn == 2,nn == 3))
        born = np.logical_and(state == 0, nn == 3)

        #Applying transitions
        state[lonely] = 0
        state[crowded] = 0
        state[alive] = 1
        state[born] = 1
        
        #Equilibrium check
        newactive = count_active(state)
        if(newactive == active):
            stepsfromeq = stepsfromeq + 1
        else:
            stepsfromeq = 0

        if((stepsfromeq == 10) & (ini == 'R')):
            eq_time = st-9
            break
        
        #'''
        #Animation
        f = open('state.dat','w')
        for i in range(r):
            for j in range(c):
                f.write('%d %d %lf\n'%(i,j,state[i,j]))
        f.close()
        plt.cla()
        im = plt.imshow(state,animated = True,vmin = 0, vmax = 1)
        plt.draw()
        plt.pause(0.001)
        #'''
        
        if(np.mod(st,50)==0):
            print('Step '+str(st)+' Completed, Active: '+str(active))
        
        st = st+1
        active = newactive

    #Converting centre of mass data
    xt = np.array(xt)
    xt =  -xt + (-1+size)
    yt = np.array(yt)
    r = np.sqrt(xt**2 + yt**2)

    if(ini == 'G'):
        fileg = open('glider.txt','w')
        fileg.write('t\tx\ty\tr\n')
        for i in range(len(comt)):
            fileg.write(str(comt[i])+'\t'+str(xt[i])+'\t'+str(yt[i])+'\t'+str(r[i])+'\n')
        fileg.close()
        
        plt.cla()
        plt.plot(comt,r)
        plt.show()

        times = comt[:90]
        rs = r[:90]
        poly = np.polyfit(times,rs,deg = 1)
    
        print('Speed: '+str(np.round(poly[0],decimals = 3))+' units per time-step')


    return eq_time#,comt,r

def hist_rand_equilibrium(n = 100,maxsteps = 5000):
    #Generate histogram of equilibriation time
    hist = []
    i = 0
    while(i<n):
        eq = main(nsteps= maxsteps,ini = 'R')
        if eq>=0:
            print('COMPLETED sweep:'+str(i))
            hist.append(eq)
            i = i+1
        else:
            continue
    

    plt.cla()
    res = plt.hist(hist,bins = 10)
    plt.show()


    vals = res[0]
    bins = res[1]
    file = open('hist.txt','w')
    file.write('value\tbin-edge\n')
    for i in range(len(vals)):
        file.write(str(vals[i])+'\t'+str(bins[i])+'\n')
    
    return hist

def com_g(xs,ys):
    #Calculate Center of Mass
    xc = np.sum(xs)/len(xs)
    yc = np.sum(ys)/len(ys)
    return xc,yc

'''
def com_glider_old(state):
    xs,ys = np.where(state == 1)
    xc = np.sum(xs)/len(xs)
    yc = np.sum(ys)/len(ys)
    return xc,yc
'''

def aliveneighbours(state):
    #Count alive neighbours
    (r,c) = state.shape
    nn = np.zeros((r,c))
    for i in range(r):
        for j in range(c):
            n = state[i,np.mod(j+1,c)]+state[np.mod(i+1,r),j]+state[i,np.mod(j-1,c)]+state[np.mod(i-1,r),j]+  state[np.mod(i+1,r),np.mod(j-1,c)]+state[np.mod(i-1,r),np.mod(j+1,c)]+state[np.mod(i-1,r),np.mod(j-1,c)]+state[np.mod(i+1,r),np.mod(j+1,c)]
            nn[i,j] = n
    return nn

def count_active(state):
    #Count total alive 
    return int(np.sum(state))

def init_state(ini,size):
    #Initialise state of board
    state = np.zeros((size,size))
    if ini == None:
        #Default
        return state
    elif ini == 'G':
        #Glider
        mid = int(np.floor(size/2))
        state[mid-1,mid-1:mid+2] = 1
        state[mid,mid+1] = 1
        state[mid+1,mid] = 1
        return state
    elif ini == 'O':
        #Oscillator
        mid = int(np.floor(size/2))
        state[mid,mid-1:mid+2] = 1
        return state
    elif ini == 'R':
        #Random
        state = np.random.randint(0,2,size = (size,size))
        return state
