import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

import matplotlib.animation as anim


def main(model = 'g',ini = '',nsteps = 10000, J = 1.0, size= 50, T = 1.0,efile = 'energies'):

    if(model not in ['g','k']):
            print('Invalid model, please set model = \'g\' for glauber and \'k\' for kawasaki')
            return
        
    acc = 0
    energies = []
    time = []
    ms = []

     
    
    spins = np.zeros((size,size))
    (r,c) = spins.shape
    if(ini ==''):
        if(model == 'g'):
            spins = np.ones((size,size))
        elif(model =='k'):
            for i in range(0,r):
                for j in range(0,c):
                    if(np.mod(i+j,2)==0):
                        spins[i,j] = +1
                    else:
                        spins[i,j] = -1
        print(spins)
    else:
        spins = ini

    #Initialise figure
    #fig = plt.figure()
    #im = plt.imshow(spins,animated=True)
        
    newspins = ''

    st = 0

    while st<=nsteps:
        if(np.mod(st,10)==0):
            print('sweep no.: '+str(st))
            #Output observables
            en = E(spins,J)
            
            energies.append(en)
            time.append(st)
            if(model =='g'):
                m = np.sum(spins)
                ms.append(m)
            #efile.write(str(st)+'\t'+str(en)+'\t'+str(m)+'\n')

            
            #Show animation
            '''
            f = open('spins.dat','w')
            for i in range(r):
                for j in range(c):
                    f.write('%d %d %lf\n'%(i,j,spins[i,j]))

            f.close()
            plt.cla()
            im = plt.imshow(spins, animated = True)
            plt.draw()
            plt.pause(0.001)
            '''
    
            
        sweep = 0
        while sweep<(size*size):
        #size*size flips each time step, each flip is random
            
            dE = 0
            f = 1
            #calc proposed change:
            if model == 'g':
                dE,newspins = glauber(spins,J)
            elif model == 'k':
                dE,newspins = kawasaki(spins,J)
                if(dE == 0):
                    sweep = sweep-1
                    f = 0
            #metropolis check
            if(f==1 & metropolis(dE,T)):
                #accept proposed change
                spins = newspins.copy()
                acc = acc+1

            sweep = sweep+1
                
        st = st+1


    efile = efile+'T'+str(np.round(T,decimals = 1))+'_model_'+model+'.dat'
    efile = open(efile,'w')
    efile.write('Time-step\tEnergy\tMagnetisation\n')   
    
    for i in range(len(time)):
        if(model =='g'):
            efile.write(str(time[i])+'\t'+str(energies[i])+'\t'+str(ms[i])+'\n')
        else:
            efile.write(str(time[i])+'\t'+str(energies[i])+'\n')
    efile.close()

    plt.cla()
    plt.plot(time, energies)
    plt.xlabel('Time-steps')
    plt.ylabel('Total Energy')
    plt.savefig('plot_T'+str(np.round(T,decimals = 1))+'_energies_model_'+model+'.png')

    print('T = '+str(T)+' completed, acc rate: '+str(np.round(100*acc/(nsteps*size*size)))+'%')
    
    return spins

def run(model,start=1.0, end=3.1, step=0.1, size = 50):
    spins = ''
    for T in np.arange(start,end, step):
        spins = main(model = model,ini = spins, T = T)

    print('Complete')
    return

def E(spins,J):
    E = 0
    (r,c) = spins.shape
    for i in range(r):
        for j in range(c):
            E = E + Eij(spins,i,j,J)/2
    return E

def Eij(spins,i,j,J):
    #energy contribution from a single spin at (i,j)
    (r,c) = spins.shape
    Eij = -J*spins[i,j]*(spins[np.mod(i+1,r),j]+spins[i,np.mod(j+1,c)]+spins[np.mod(i-1,r),j]+spins[i,np.mod(j-1,c)])

    return Eij

def dE_G(spins,newspins,i,j,J):
    #change in energy in glauber dynamics
    
    (r,c) = spins.shape
    Ef = Eij(newspins,i,j,J)
    Ei = Eij(spins,i,j,J)
    
    return Ef-Ei

def glauber(spins,J):
    #glauber dynamics: flip spin
    
    dE = 0
    (r,c) = spins.shape

    #select random index
    pos = np.random.random(2)*np.array([r,c])
    (i,j) = int(np.floor(pos[0])),int(np.floor(pos[1]))

    #copy to a new array
    newspins = spins.copy()

    #flip the spin
    newspins[i,j] = -newspins[i,j]

    #calc change in energy
    dE = dE_G(spins,newspins,i,j,J)
    
    return dE,newspins

def dE_K(spins,newspins,i1,j1,i2,j2,J):
    #change in energy in kawasaki dynamics
    
    (r,c) = spins.shape
    Ef = Eij(newspins,i1,j1,J)+Eij(newspins,i2,j2,J)
    Ei = Eij(spins,i1,j1,J)+Eij(spins,i2,j2,J)
    
    return Ef-Ei

def kawasaki(spins,J):
    #kawasaki dynamics: swap two spins
    dE = 0
    (r,c) = spins.shape

    #select random index
    pos = np.random.random(4) *np.array([r,c,r,c])
    (i1,j1,i2,j2) = int(np.floor(pos[0])),int(np.floor(pos[1])),int(np.floor(pos[2])),int(np.floor(pos[3]))

    if((i1==i2 & j1==j2)|(spins[i1,j1] == spins[i2,j2])):
        #checking if change irrelevant(equal spins/same index)
        return dE, spins
    else:
        #swapping the spins
        newspins = spins.copy()
        newspins[i1,j1] = spins[i2,j2]
        newspins[i2,j2] = spins[i1,j1]

        #calc change in energy from swap
        dE = dE_K(spins,newspins,i1,j1,i2,j2,J)

        return dE,newspins

def metropolis(dE,T):
    #metropolis condition to accept or reject proposed change
    
    flag = False
    kB = 1.0
    metfac = np.exp(-dE/(kB*T))
    if(metfac<1):
        ran = np.random.random()
        if(ran<=metfac):
            flag = True
    else:
        flag = True
    return flag

