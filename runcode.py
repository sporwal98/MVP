from mcising import *

def run(model,start=1.0, end=3.1, step=0.1):
    for T in np.arange(start,end, step):
        main(model = model, T = T)

    print('Complete')

    return
