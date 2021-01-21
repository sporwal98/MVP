import numpy as np
import sys
import matplotlib.pyplot as plt
import sys

def main(model = 'g',nsteps = 10000, J = 1.0, size = 50, T = 1.0,efile = 'energies'):
    acc = 0
    energies = []
    time = []

    efile = efile+'T'+str(T)+'_model_'+model+'.txt'
    efile = open(efile,'w')
    efile.write('Time-step\tEnergy\n')
    
    
    spins = np.zeros((size,size))
    (r,c) = spins.shape
    #randomly initialise spins
    for i in range(r):
        for j in range(c):
            ran = np.random.random()
            if(ran>0.5):
                spins[i,j] = +1
            else:
                spins[i,j] = -1
    print(spins)
    
    newspins = spins.copy()

    for st in range(nsteps):
        if(np.mod(st,10)==0):
            print('sweep no.: '+str(st))
            en = E(spins,J)
            energies.append(en)
            time.append(st)
            efile.write(str(st)+'\t'+str(en)+'\n')
        
        for sweep in range(size*size):
        #size*size flips each time step, each flip is random
            dE = 0
            #calc proposed change:
            if model == 'g':
                dE,newspins = glauber(spins,J)
            else:
                dE,newspins = kawasaki(spins,J)

            #metropolis check
            if(metropolis(dE,T)):
                spins = newspins.copy()
                acc = acc+1

    efile.close()
    plt.cla()
    plt.plot(time, energies)
    plt.xlabel('Time-steps')
    plt.ylabel('Total Energy')
    
    plt.savefig('plot_T'+str(T)+'_energies.png')
    print('T = '+str(T)+' completed')   
    #print('acc rate: '+str(100*acc/(nsteps*size*size)))
    
    return

def E(spins,J):
    E = 0
    (r,c) = spins.shape
    for i in range(r):
        for j in range(c):
            E = E + Eij(spins,i,j,J)
    return E

def Eij(spins,i,j,J):
    (r,c) = spins.shape
    Eij = -J*spins[i,j]*(spins[np.mod(i+1,r),j]+spins[i,np.mod(j+1,c)]+spins[np.mod(i-1,r),j]+spins[i,np.mod(j-1,c)])/2

    return Eij

def dE_G(spins,newspins,i,j,J):
    (r,c) = spins.shape
    Ef = Eij(newspins,i,j,J)
    Ei = Eij(spins,i,j,J)
    
    return Ef-Ei

def glauber(spins,J):
    #flip spin
    dE = 0
    (r,c) = spins.shape
    pos = np.random.random(2) *np.array([r-1,c-1])

    (i,j) = int(np.round(pos[0])),int(np.round(pos[1]))
    newspins = spins.copy()
    newspins[i,j] = -newspins[i,j]

    dE = dE_G(spins,newspins,i,j,J)
    return dE,newspins

def dE_K(spins,newspins,i1,j1,i2,j2,J):
    (r,c) = spins.shape
    Ef = Eij(newspins,i1,j1,J)+Eij(newspins,i2,j2,J)
    Ei = Eij(spins,i1,j1,J)+Eij(spins,i2,j2,J)
    
    return Ef-Ei

def kawasaki(spins,J):
    #swap spins
    dE = 0
    (r,c) = spins.shape
    pos = np.random.random(4) *np.array([r-1,c-1,r-1,c-1])

    (i1,j1,i2,j2) = int(np.round(pos[0])),int(np.round(pos[1])),int(np.round(pos[2])),int(np.round(pos[3]))
    if((i1==i2 & j1==j2)|(spins[i1,j1] == spins[i2,j2])):
        newspins = spins.copy()
    else:
        newspins = spins.copy()
        newspins[i1,j1] = spins[i2,j2]
        newspins[i2,j2] = spins[i1,j1]
        dE = dE_K(spins,newspins,i1,j1,i2,j2,J)

    return dE,newspins

def metropolis(dE,T):
    flag = False
    kB = 1.0
    metfac = np.exp(-dE/(kB*T))
    if(metfac<1):
        ran = np.random.random()
        if(ran<metfac):
            flag = True
    else:
        flag = True
    return flag

