import matplotlib
matplotlib.use('TKAgg')

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial

import matplotlib.animation as anim
import pandas as pd

const_k = 1/(4*np.pi)#epsilon0 = 0
def SOR(wmin = 1.0, wmax = 2.0,dw = 0.01):
    ws = np.arange(wmin, wmax, dw)
    itr = []
    for w in ws:
        print('w:',w)
        itr.append(main(w = w, update = 'SOR', threshold = 1e-2))
        
    df = pd.DataFrame.from_dict({'w':ws,'iterations':itr})
    df.to_csv('itrvsws_'+str(wmin)+'_'+str(wmax)+'.csv')

    return ws, itr

def main(w = 1.94,size = 50,dx = 1, dist = 'PT',update = 'J',threshold = 1e-3, file = 'data'):
    #Code to run simulation
    
    rho = np.zeros((size,size,size))
    if(dist == 'PT'):
        rho[int(size/2),int(size/2),int(size/2)] = 1
    elif(dist == 'WR'):
        rho[int(size/2),int(size/2),:] = 1

    phi = rho/6

    st = 0
    while(True):
        fixboundary(phi)
        if(update == 'J'):
            #Jacobi Update
            phinew = np.copy(phi)
            phinew= (np.roll(phi,1,axis = 0)+np.roll(phi,-1,axis = 0)+
                      np.roll(phi,1,axis = 1)+np.roll(phi,-1,axis = 1)+
                      np.roll(phi,1,axis = 2)+np.roll(phi,-1,axis = 2)+
                      rho)/6
            fixboundary(phinew)
            diff = np.sum(np.abs(phinew-phi))
            phi = np.copy(phinew)

        elif(update == 'GS'):
            #Gauss- Seidel
            diff = 0
            for i in range(1,size-1):
                for j in range(1,size-1):
                    for k in range(1,size-1):
                        tmp = phi[i,j,k]
                        phi[i,j,k] = (phi[i+1,j,k]+phi[i-1,j,k]+
                                      phi[i,j+1,k]+phi[i,j-1,k]+
                                      phi[i,j,k+1]+phi[i,j,k-1]+
                                      rho[i][j][k])/6
                        diff += np.abs(phi[i,j,k]-tmp)

        elif(update == 'SOR'):
            #GS w SOR
            diff = 0
            for i in range(1,size-1):
                for j in range(1,size-1):
                    for k in range(1,size-1):
                        tmp = phi[i,j,k]
                        phi[i,j,k] = w*(phi[i+1,j,k]+phi[i-1,j,k]+
                                      phi[i,j+1,k]+phi[i,j-1,k]+
                                      phi[i,j,k+1]+phi[i,j,k-1]+
                                      rho[i][j][k])/6 + (1-w)*phi[i,j,k]
                        diff += np.abs(phi[i,j,k]-tmp)

        if(np.mod(st,100) == 0):
            print(st,', diff ' ,diff)
        if(diff<threshold):
            break
        st +=1

        
    '''
    plt.cla()
    plt.contourf(phi[:,:,int(size/2)])
    plt.colorbar()
    plt.title('Countour plot of potential along cross-section('+dist+')')
    plt.savefig('Potential_Cross-section_'+update+'_dist_'+dist+'.png')
    #'''

    '''
    Ex,Ey,Ez = (calcEx(phi,dx),calcEy(phi,dx),calcEz(phi,dx))
    E = np.sqrt(Ex**2+Ey**2+Ez**2)
    normEx,normEy,normEz = Ex/E,Ey/E,Ez/E
    #'''

    '''
    plt.figure()
    plt.cla()
    plt.quiver(normEy[:,:,int(size/2)],normEx[:,:,int(size/2)])
    plt.savefig('E-field_'+update+'dist_'+dist+'.png')
    #'''
    
    '''
    (Bx,By,Bz) = calcB(phi,dx)
    B = np.sqrt(Bx**2+By**2+Bz**2)
    normBx,normBy,normBz = Bx/B,By/B,Bz/B 
    #'''
    '''
    plt.figure()
    plt.cla()
    plt.quiver(normBy[:,:,int(size/2)],normBx[:,:,int(size/2)])
    plt.savefig('B-field_'+update+'_dist_'+dist+'.png')
    #'''

    #'''
    rowi = []
    colj = []
    depk = []
    rs = []
    pots = []
    Exs = []
    Eys = []
    Ezs = []
    Bxs = []
    Bys = []
    Bzs = []
    #'''

    '''
    #3D
    for i in range(1,size-1):
        for j in range(1,size-1):
            for k in range(1,size-1):
                r = np.sqrt((i-(size/2))**2+(j-(size/2))**2+(k-(size/2))**2)
                #efile.write(str(i)+'\t'+str(j)+'\t'+str(k)+'\t'+str(r)+'\t'
                #           str(phi[i,j,k])+'\t'+
                #           str(Ex[i,j,k])+'\t'+str(Ey[i,j,k])+'\t'+str(Ez[i,j,k])+'\n')

                rowi.append(i)
                colj.append(j)
                depk.append(k)
                rs.append(r)
                pots.append(phi[i,j,k])
                Exs.append(Ex[i,j,k])
                Eys.append(Ey[i,j,k])
                Ezs.append(Ez[i,j,k])
                Bxs.append(Bx[i,j,k])
                Bys.append(By[i,j,k])
                Bzs.append(Bz[i,j,k])
                
    df = pd.DataFrame.from_dict({
            'i':rowi,
            'j':colj,
            'k':depk,
            'r':rs,
            'Potential':pots,
            'Ex':Exs,
            'Ey':Eys,
            'Ez':Ezs,
            'Bx':Bxs,
            'By':Bys,
            'Bz':Bzs
        })
    df.to_csv('file_update_'+update+'_tol_'+str(threshold)+'_dist_'+dist+'.csv')
    #'''

    #2D along cut
    '''
    rowi = []
    colj = []
    depk = []
    rs = []
    pots = []
    Exs = []
    Eys = []
    Ezs = []
    Bxs = []
    Bys = []
    Bzs = []
    for i in range(1,size-1):
        for j in range(1,size-1):
            k = int(size/2)
            r = np.sqrt((i-(size/2))**2+(j-(size/2))**2+(k-(size/2))**2)
            #efile.write(str(i)+'\t'+str(j)+'\t'+str(k)+'\t'+str(r)+'\t'
            #           str(phi[i,j,k])+'\t'+
            #           str(Ex[i,j,k])+'\t'+str(Ey[i,j,k])+'\t'+str(Ez[i,j,k])+'\n')

            rowi.append(i)
            colj.append(j)
            depk.append(k)
            rs.append(r)
            pots.append(phi[i,j,k])
            Exs.append(Ex[i,j,k])
            Eys.append(Ey[i,j,k])
            Ezs.append(Ez[i,j,k])
            Bxs.append(Bx[i,j,k])
            Bys.append(By[i,j,k])
            Bzs.append(Bz[i,j,k])
                
    df = pd.DataFrame.from_dict({
            'i':rowi,
            'j':colj,
            'k':depk,
            'r':rs,
            'Potential':pots,
            'Ex':Exs,
            'Ey':Eys,
            'Ez':Ezs,
            'Bx':Bxs,
            'By':Bys,
            'Bz':Bzs
        })
    df.to_csv('file_update_'+update+'_tol_'+str(threshold)+'_dist_'+dist+'2D.csv')
    

    
    #efile.close()
    #'''
    
    
    print('concluded')
    return st


def behaviourwithr(file):
    from numpy.polynomial import polynomial as poly
    df = pd.read_csv(file)
    r = df.r.to_numpy()
    U = df.Potential.to_numpy()
    E = np.sqrt((df.Ex.to_numpy())**2+(df.Ey.to_numpy())**2+(df.Ez.to_numpy())**2)
    B = np.sqrt((df.Bx.to_numpy())**2+(df.By.to_numpy())**2+(df.Bz.to_numpy())**2)

    filt = r!=0
    r = r[filt]
    U = U[filt]
    E = E[filt]
    B = B[filt]

    
    plt.cla()
    plt.scatter(r,U,marker = 'x')
    plt.xlabel('r')
    plt.ylabel('U')
    plt.savefig('UvR'+file+'.png')

    plt.cla()
    plt.scatter(r,E,marker = 'x')
    plt.xlabel('r')
    plt.ylabel('E')
    plt.savefig('EvR'+file+'.png')

    plt.cla()
    plt.scatter(r,B,marker = 'x')
    plt.xlabel('r')
    plt.ylabel('B')
    plt.savefig('BvR'+file+'.png')


    
    #POINT CHARGE
    #FIT TO RANGE 0.5-1.5
    fitx = np.log(r[U!=0])
    fity = np.log(U[U!=0])
    fitrange = np.logical_and(fitx>0.5,fitx<1.5)
    fitx = fitx[fitrange]
    fity = fity[fitrange]
    Ufit = poly.polyfit(fitx,fity,deg = 1)
    Upred = Ufit[0] + np.log(r[U!=0])*Ufit[1]
    plt.cla()
    plt.scatter(np.log(r[U!=0]),np.log(U[U!=0]),marker = 'x')
    plt.plot(np.log(r[U!=0]),Upred,label = 'FIT: slope = '+str(Ufit[1]),c = 'r')
    plt.xlabel('logR')
    plt.ylabel('logU')
    plt.legend()
    plt.savefig('logUvlogR'+file+'.png')

    #FIT RANGE 0.5-2.5
    fitx = np.log(r[E!=0])
    fity = np.log(E[E!=0])
    fitrange = np.logical_and(fitx>0.5,fitx<2.5)
    fitx = fitx[fitrange]
    fity = fity[fitrange]
    Efit = poly.polyfit(fitx,fity,deg = 1)
    Epred = Efit[0] + np.log(r[E!=0])*Efit[1]
    plt.cla()
    plt.scatter(np.log(r[E!=0]),np.log(E[E!=0]),marker = 'x')
    plt.plot(np.log(r[E!=0]),Epred,label = 'FIT: slope = '+str(Efit[1]),c = 'r')
    plt.xlabel('logR')
    plt.ylabel('logE')
    plt.legend()
    plt.savefig('logEvlogR'+file+'.png')



    #WIRE
    plt.cla()
    plt.scatter(np.log(r[U!=0]),U[U!=0],marker ='x')
    plt.xlabel('logR')
    plt.ylabel('U')
    plt.legend()
    plt.savefig('UvlogR'+file+'.png')
    
    #FIT RANGE 1-3
    fitx = np.log(r[B!=0])
    fity = np.log(B[B!=0])
    fitrange = np.logical_and(fitx>1,fitx<3)
    fitx = fitx[fitrange]
    fity = fity[fitrange]
    Bfit = poly.polyfit(fitx,fity,deg = 1)
    Bpred = Bfit[0] + np.log(r[B!=0])*Bfit[1]
    plt.cla()
    plt.scatter(r[B!=0],np.log(B[B!=0]),marker ='x')
    plt.plot(np.log(r[B!=0]),Bpred,label = 'FIT: slope = '+str(Bfit[1]),c = 'r')
    plt.legend()
    plt.xlabel('logR')
    plt.ylabel('logB')
    plt.savefig('logBvlogR'+file+'.png')

def itrvswplot(file):
    df = pd.read_csv(file)
    ws = df.w.to_numpy()
    itr = df.iterations.to_numpy()

    plt.cla()
    plt.scatter(ws,itr,marker = 'x')
    plt.xlabel('w')
    plt.ylabel('iterations')
    plt.savefig('itrvsws.png')
    

def fixboundary(a):
    #set boundary to 0
    for i in range(-1,1):
        for j in range(-1,1):
            for k in range(-1,1):
                a[i,:,:]= 0
                a[:,j,:]= 0
                a[:,:,k] = 0


    
def calcEx(phi,dx):
    #Ex = -dphi/dx
    return -(np.roll(phi,-1,axis = 0)-np.roll(phi,+1,axis = 0))/(2*dx)

def calcEy(phi,dx):
    #Ey = -dphi/dy
    return -(np.roll(phi,-1,axis = 1)-np.roll(phi,+1,axis = 1))/(2*dx)

def calcEz(phi,dx):
    #Ez = -dphi/dz
    return -(np.roll(phi,-1,axis = 2)-np.roll(phi,+1,axis = 2))/(2*dx)

def calcB(phi,dx):
    #B = curl(vec(V) =  {0,0,V})
    Bx = (np.roll(phi,1,axis = 1)-np.roll(phi,-1,axis = 1))/(2*dx)
    By = -(np.roll(phi,1,axis = 0)-np.roll(phi,-1,axis = 0))/(2*dx)
    Bz = np.zeros(Bx.shape)

    return (Bx,By,Bz)
