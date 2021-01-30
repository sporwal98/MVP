import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main(model):
    Ts = np.round(np.arange(1.0,3.1,0.1),decimals = 1)#T = 1.0 NEEDS CHECKING for model G
    dfs = []
    avgens = []
    avgms = []
    sphs = []
    sus = []
    
    

    
    for i in range(len(Ts)):
        fname = 'energiesT'+str(Ts[i])+'_model_'+model+'.dat'
        df = pd.read_table(fname)
        dfs.append(df)

        avgens.append(np.mean( (df['Energy'])[100:] ))
        sphs.append(specific_heat((df['Energy'])[100:],Ts[i]))
        if(model =='g'):
            avgms.append(np.mean((df['Magnetisation'])[100:]))
            sus.append(susceptibility((df['Magnetisation'])[100:],Ts[i]))
        

    print(sus[:10])
    print(sphs[:10])

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
    plt.plot(Ts,sphs)
    plt.xlabel('T')
    plt.ylabel('Specific Heat')
    plt.title('Specific Heat vs T')
    #plt.show()
    plt.savefig('CvT_'+model+'.png')
    
    file = 'observables_'+model+'.dat'
    file = open(file,'w')
    file.write('E\tC\tM\tS\n')
    for i in range(len(Ts)):
        if model =='g':
            file.write(str(avgens[i])+'\t'+
                       str(sphs[i])+'\t'+
                       str(avgms[i])+'\t'+
                       str(sus[i])+'\n')
        else:
            file.write(str(avgens[i])+'\t'+
                       str(sphs[i])+'\n')

    file.close()
    return

def susceptibility(mags,T,N=2500, kB = 1.0):
    sus = (np.mean(mags**2)-np.mean(mags)**2)/(N*kB*T)
    return sus

def specific_heat(ens,T,N = 2500, kB = 1.0):
    sph = (np.mean(ens**2)-np.mean(ens)**2)/(N*kB*T**2)
    return sph
