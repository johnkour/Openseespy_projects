"""
A basic program to plot in the same PDF file P - delta curves from text files 
generated from Openseespy recorders.

@author: John Kouretas
"""


##############################################################################


import numpy as np
import matplotlib.pyplot as plt


##############################################################################


def Plot_Pdelta(F, filename):
    
    fig, ax = plt.subplots()
    
    ax.set_xlabel('Vertical displacement at the middle [dm]')
    ax.set_ylabel('Force [MN]')
    
    ax.set_title('P-delta curves for vertical imperfection equal to 10% x Pcr')
    plt.grid(linestyle='--', linewidth=2)
    
    for Analysis in ['Linear', 'PDelta', 'Corotational']:
        
        fname = Analysis + '/'
        fname += filename
        temp = np.loadtxt(fname, delimiter = ' ')
        Data = np.column_stack((temp[:, 0], temp[:, 2]))
        #Data = np.insert(Data, [0], [[0, 0]], axis= 0)
        
        P = Data[:, 0]
        P *= F
        P /= 10**6
            
        disp = Data[:, 1]
        disp *= 1
        disp *= 10**1
        
        lnstyle = 'solid'
        if Analysis == 'PDelta':
            lnstyle = 'dotted'
        #(0, (5, 10))
        
        ax.plot(disp, P, label = Analysis + ' ' + 'analysis', 
                linestyle = lnstyle)
        
    Data = np.loadtxt('EulerLoad.txt', delimiter = ' ')
    
    
    P = Data[:, 0]
    P /= 10**6
        
    disp = Data[:, 1]
    disp *= 10**1
    
    lnstyle = 'dashdot'
    
    ax.plot(disp, P, label = "Euler's Buckling Load", linestyle = lnstyle)
    
    #ax.legend(loc = 'lower left')
    ax.legend(loc = 'upper right')
    
    ax.set_xlim(left = 0, right = 55)
    #ax.set_ylim(bottom = 0, top = 100)
    
    filename = filename[:-4]
    filename += '.png'
    plt.savefig(filename)
    
    filename = filename[:-4]
    filename += '.pdf'
    plt.savefig(filename)
    
    plt.show()
    

##############################################################################


text1 = '''
Please enter the name of the file containing the results of the analysis (don't write the ".txt"):
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


