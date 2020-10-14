"""
A basic program to plot P - delta curves from text files generated from
Openseespy recorders.

@author: John Kouretas
"""


##############################################################################


import numpy as np
import matplotlib.pyplot as plt


##############################################################################


def Plot_Pdelta(F, filename):
    
    
    
    Data = np.loadtxt(filename, delimiter = ' ')
    np.insert(Data, [0], [[0, 0]], axis= 0)
    
    fig, ax = plt.subplots()
    
    P = Data[:, 0]
    P *= -F
    P /= 10**3
    
    disp = Data[:, 1]
    disp *= -1
    disp *= 100
    
    ax.plot(disp, P)
    
    ax.set_xlabel('displacement [cm]')
    ax.set_ylabel('Force [kN]')
    
    ax.set_title('P-delta curve')
    plt.grid(linestyle='--', linewidth=2)
    
    filename = filename[:-4]
    filename += '.png'
    plt.savefig(filename)
    
    filename = filename[:-4]
    filename += '.pdf'
    plt.savefig(filename)
    
    plt.show()
    

##############################################################################


text1 = '''
Please enter the path to the file you want to plot (don't write the ".txt"):
'''

filename = input(text1)
temp = filename
filename += '.txt'

text2 = '''
Please enter the value of the load that was applied when analizing in N (P):
'''

P = input(text2)
P = float(P)

Plot_Pdelta(P, filename)


##############################################################################


