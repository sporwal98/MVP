import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main(model):
    #Generates the Temperature plots for Avg Energy, Avg Magnetisation(if applicable), Specific Heat and Susceptibility(if applicable)
    '''
    Inputs:
    model = 'g' for Glauber and 'k' for Kawasaki
    '''
    Ts = np.round(np.arange(1.0,3.1,0.1),decimals = 1)#T = 1.0 NEEDS CHECKING for model G
    dfs = []
    avgens = []
    avgms = []
    sphs = []
    sus = []
    
    avgenerr = []
    for i in range(len(Ts)):
        #Reading energies and Magnetisation from data files, generating errors
        fname = 'energiesT'+str(Ts[i])+'_model_'+model+'.dat'
        df = pd.read_table(fname)
        dfs.append(df)

        avgens.append(np.mean( (df['Energy'])[100:] ))
        sphs.append(specific_heat((df['Energy'])[100:],Ts[i]))
        if(model =='g'):
            avgms.append(np.mean((df['Magnetisation'])[100:]))
            sus.append(susceptibility((df['Magnetisation'])[100:],Ts[i]))
            
        avgenerr.append(jackknife(df['Energy'][100:]))
        print('Read Complete:' + str(Ts[i]))

    print(sus[:10])
    print(sphs[:10])

    avgens = np.array(avgens)
    avgms = np.array(avgms)
    sphs = np.array(sphs)
    sus = np.array(sus)
    avgenerr = np.array(avgenerr)

    plt.cla()
    plt.plot(Ts,avgens)
    plt.xlabel('T')
    plt.ylabel('Average energy')
    plt.title('Avg Energy vs T')
    #plt.show()
    plt.savefig('EvT_'+model+'.png')

    if(model == 'g'):
        plt.cla()
        plt.plot(Ts,np.abs(avgms))
        plt.xlabel('T')
        plt.ylabel('Average Magnetisation')
        plt.title('Avg Magnetisation vs T')
        #plt.show()
        plt.savefig('MvT_'+model+'.png')
    
        plt.cla()
        plt.plot(Ts,sus)
        plt.xlabel('T')
        plt.ylabel('Susceptibility')
        plt.title('Susceptibility vs T')
        #plt.show()
        plt.savefig('SvT_'+model+'.png')

    plt.cla()
    plt.plot(Ts,sphs, label = 'Value')
    plt.plot(Ts,sphs+avgenerr,label = 'Positive Error')
    plt.plot(Ts,sphs-avgenerr, label = 'Negative Error')
    plt.xlabel('T')
    plt.ylabel('Specific Heat')
    plt.title('Specific Heat vs T')
    #plt.show()
    plt.savefig('CvT_'+model+'.png')
    
    file = 'observables_'+model+'.dat'
    file = open(file,'w')
    file.write('T\tE\tC\tdeltaC\tM\tS\n')
    for i in range(len(Ts)):
        if model =='g':
            file.write(str(Ts[i])+'\t'+
                       str(avgens[i])+'\t'+
                       str(sphs[i])+'\t'+
                       str(avgenerr[i])+'\t'+
                       str(avgms[i])+'\t'+
                       str(sus[i])+'\n')
        else:
            file.write(str(Ts[i])+'\t'+
                       str(avgens[i])+'\t'+
                       str(sphs[i])+'\t'+
                       str(avgenerr[i])+'\n')

    file.close()
    return

def susceptibility(mags,T,N=2500, kB = 1.0):
    #Calculate susceptibility given magnetisations
    sus = (np.mean(mags**2)-np.mean(mags)**2)/(N*kB*T)
    return sus

def specific_heat(ens,T,N = 2500, kB = 1.0):
    #Calculate specific heat given energies
    sph = (np.mean(ens**2)-np.mean(ens)**2)/(N*kB*T**2)
    return sph

def jackknife(data):
    #Calculate error using jackknife resampling
    #print(data[:10])
    n = len(data)
    xibar = []
    for i in range(n):
        xi = 0
        for j in range(n):
            if(i!=j):
                xi = xi + data[100+j]
        xi = xi/(n-1)
        xibar.append(xi)
    xibar = np.array(xibar)
    
    mean = np.mean(xibar)

    var = np.sum((xibar-mean)**2)*(n-1)/n

    #print(var)
    return np.sqrt(var)/np.sqrt(n)
