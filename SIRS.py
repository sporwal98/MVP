import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

import matplotlib.animation as anim
def main(p1,p2,p3,imm = 0,size= 50, nsteps = 1100):
    #SIRS Simulation
    states = ['S','I','R']
    
    grid = grid_init(size,imm)

    print('S: '+str(np.sum(grid == 0))+',I: '+str(np.sum(grid == 1))+',R: '+str(np.sum(grid ==2))+',Im: '+str(np.sum(grid == 3))+' ('+str(imm)+')')
    r,c = grid.shape
    N = r*c
    
    #'''
    #Figure Setups
    fig = plt.figure()
    im = plt.imshow(grid, animated = True,vmin = 0,vmax = 3)
    plt.colorbar(im)
    #'''
    
    infected = []

    st = 0
    while(st<nsteps):
        #One major sweep
        if(st>99):
            print(count_I(grid))
            #After 100 equilibriation steps
            infected.append(count_I(grid))
            #print(infected)
            #'''
            if(infected[-1] == 0):
                #Checking cured
                while(st<nsteps):
                    infected.append(0)
                    st = st+1
                return 0.0,0.0
            #'''

        #'''
        #Animation
        f = open('grid.dat','w')
        for i in range(r):
            for j in range(c):
                f.write('%d %d %lf\n'%(i,j,grid[i,j]))
        f.close()
        plt.cla()
        im = plt.imshow(grid,animated = True,vmin = 0, vmax = 3)
        plt.draw()
        plt.pause(0.001)
        #'''
        

        swst = 0
        while(swst<N):
            #Each step in the sweep
            i = np.random.randint(0,r)
            j = np.random.randint(0,c)
            randcheck = np.random.random()
            
            if((grid[i,j] == 0) & infected_neighbour(grid,i,j) & (randcheck<p1)):
                grid[i,j] = 1
            elif((grid[i,j] == 1) & (randcheck<p2)):
                grid[i,j] = 2
            elif((grid[i,j] == 2) & (randcheck<p3)):
                grid[i,j] = 0
            
            swst = swst + 1
        if(np.mod(st,500) == 0):
            print('Sweep no. '+str(st)+' Completed')
        st = st+1

    

    #print(infected)

    
    return (np.mean(infected)/N), (np.var(infected)/N)

    
        
        

def grid_init(size,imm):
    #Initialise grid
    grid = np.zeros((size,size))
    if(imm == 0):
        #init without immunity
        grid = np.random.randint(0,3,(size,size))
    else:
        #init with immunity
        r,c = grid.shape
        for i in range(r):
            for j in range(c):
                rand = np.random.random()
                if(rand<imm):
                    grid[i][j] = 3
                else:
                    grid[i][j] = np.random.randint(0,3)

    return grid

def infected_neighbour(grid,i,j):
    #Infected neighbour check
    r,c = grid.shape
    return (1 in [grid[np.mod(i+1,r),j],grid[np.mod(i-1,r),j],grid[i,np.mod(j+1,c)],grid[i,np.mod(j-1,c)]])

def count_I(grid):
    #Number of infected
    return np.sum(grid == 1)

def phasediagram_3(p2,reso=0.05):
    #Generate phase diagram data p2 = 0.5
    vals = np.round(np.arange(0,1.00001,reso),decimals = 2)
    file = open('phasediag.txt',mode = 'w')
    file.write('p1\tp3\t<I>/N\tVar(I)\n')
    
    
    cont = np.zeros((len(vals),len(vals)))
    for i in range(len(vals)):
        for j in range(len(vals)):
            p1 = vals[i]
            p3 = vals[j]
            print('p1: '+str(p1)+',p3: '+str(p3))
            meanI,varI = main(p1,p2,p3,size = 50,nsteps = 1100)
            cont[i,j] = meanI
            file.write(str(p1)+'\t'+str(p3)+'\t'+str(meanI)+'\t'+str(varI)+'\n')
    file.close()
    '''
    plt.cla()
    fig = plt.figure()
    im = plt.imshow(cont,vmin=0,vmax = 1)
    
    plt.colorbar(im)
    plt.draw()
    plt.show()
    #'''
    
def contourplotfromfile_mean(fname,size = 50):
    #phasediag.txt plot mean
    file = open(fname, mode = 'r')

    head = file.readline()
    p1s = []
    p3s = []
    meanIs = []

    while(True):
        line = file.readline()
        line = line.split('\t')
        if(line[0]==''):
            break
        line[3] = (line[3])[:-1]
        #print(line)
        p1s.append(line[0])
        p3s.append(line[1])
        meanIs.append(line[2])


    print(meanIs)
    p1s = np.array(p1s)
    p3s = np.array(p3s)
    meanIs = np.array(meanIs)
    
    vals = np.round(np.arange(0,1.00001,0.05),decimals = 2)
    cont = np.zeros((len(vals),len(vals)))
    m = 0
    for i in range(len(vals)):
        for j in range(len(vals)):
            cont[i,j] = meanIs[m]
            m+=1
    #print(cont)
    plt.cla()
    plt.contourf(vals,vals,cont)
    #plt.imshow(cont)
    plt.xlabel('p1')
    plt.ylabel('p3')
    plt.colorbar()

    plt.show()

def contourplotfromfile_var(fname,size = 50):
    #phasediag.txt plot var
    file = open(fname, mode = 'r')

    head = file.readline()
    p1s = []
    p3s = []
    meanIs = []
    varIs = []

    while(True):
        line = file.readline()
        line = line.split('\t')
        if(line[0]==''):
            break
        line[3] = (line[3])[:-1]
        #print(line)
        p1s.append(line[0])
        p3s.append(line[1])
        meanIs.append(line[2])
        varIs.append(line[3])

    p1s = np.array(p1s)
    p3s = np.array(p3s)
    meanIs = np.array(meanIs)
    varIs = np.array(varIs)
    
    vals = np.round(np.arange(0,1.00001,0.05),decimals = 2)
    cont = np.zeros((len(vals),len(vals)))
    m = 0
    for i in range(len(vals)):
        for j in range(len(vals)):
            cont[i,j] = varIs[m]
            m+=1
    #print(cont)
    plt.cla()
    plt.contourf(vals,vals,cont)
    #plt.imshow(cont)
    plt.xlabel('p1')
    plt.ylabel('p3')
    plt.colorbar()

    plt.show()

def variancealongfixedcut(fname = 'savefixedcut.txt',size= 50):
    #Finding variance in p2 = 0.5,p3 = 0.5, p1 in range 0.2 to 0.5
    mea = []
    var = []
    p2 = 0.5
    p3 = 0.5
    p1s = np.round(np.arange(0.2,0.55,0.05),decimals = 2)
    for p1 in p1s:
        #print('p1: '+str(p1)+'p2: '+str(p2)+',p3: '+str(p3))
        mI, vI = main(p1,p2,p3, size = size, nsteps = 10100)
        mea.append(mI)
        var.append(vI)

    file = open('fixedcut_p2_'+str(p2)+'_p3_'+str(p3)+'.txt','w')
    file.write('p1\tp2\tp3\tmeanI\tmeanIerr\tvarI\tvarIerr\n')
    meanIerr = jackknife(mea)
    varIerr = jackknife(var)
    for i in range(len(p1s)):
        file.write(str(p1s[i])+'\t'+str(p2)+'\t'+str(p3)+'\t'+str(mea[i])+'\t'+str(meanIerr)+'\t'+str(var[i])+'\t'+str(varIerr)+'\n')

    var = np.array(var)
    
    file.close()
    plt.cla()
    plt.plot(p1s,var)
    plt.plot(p1s,var+varIerr)
    plt.plot(p1s,var-varIerr)
    plt.xlabel('p1')
    plt.ylabel('variance')
    plt.title('p2 = '+str(p2)+',p3 = '+str(p3))
    plt.savefig('fixedcutplot.png')

    return p1s,mea,var

def fracimmune(p1,p2,p3,size = 50):
    #Immune fraction run
    imm = np.arange(0.0,1.001,0.01)
    brk = 0
    meanIs = []
    varIs  = []
    for i in range(len(imm)):
        mI, vI = main(p1,p2,p3,imm = imm[i], size = size)
        
        meanIs.append(mI)
        varIs.append(vI)

        if(mI == 0):
            brk = imm[i]
            break
    for i in range(len(meanIs),len(imm)):
        meanIs.append(0.0)
        varIs.append(0.0)
        
    meanIerr = jackknife(meanIs)
    varIerr = jackknife(varIs)
    file = 'fracImmunePlot_'+str(p1)+'_'+str(p2)+'_'+str(p3)+'_'+str(size)+'.txt'
    file = open(file,'w')
    file.write('imm\tmeanI\tmeanIerr\tvarI\tvarIerr\n')
    for i in range(len(meanIs)):
        file.write(str(imm[i])+'\t'+str(meanIs[i])+'\t'+str(meanIerr)+'\t'+str(varIs[i])+'\t'+str(varIerr)+'\n')
    file.close()   

    meanIs = np.array(meanIs)
    meanIerr = np.array(meanIerr)
    plt.cla()
    plt.plot(imm,meanIs)
    plt.plot(imm,meanIs+meanIerr)
    plt.plot(imm,meanIs-meanIerr)

    plt.xlabel('immune ratio')
    plt.ylabel('mean I')
    plt.savefig('Immratio_plot_'+str(p1)+'_'+str(p2)+'_'+str(p3)+'cured_at_imm_'+str(brk)+'.png')
    
    
    return

def jackknife(data):
    #Calculate error using jackknife resampling
    
    n = len(data)
    xibar = []
    for i in range(n):
        xi = 0
        for j in range(n):
            if(i!=j):
                xi = xi + data[j]
        xi = xi/(n-1)
        xibar.append(xi)
    xibar = np.array(xibar)
    
    mean = np.mean(xibar)

    var = np.sum((xibar-mean)**2)*(n-1)/n

    #print(var)
    return np.sqrt(var)/np.sqrt(n)
        
    
    
    
    
    
    
