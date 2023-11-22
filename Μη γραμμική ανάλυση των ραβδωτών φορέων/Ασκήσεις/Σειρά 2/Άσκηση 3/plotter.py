"""
A basic program to plot in the same PDF file P - delta curves from text files 
generated from Openseespy recorders.

@author: John Kouretas
"""


##############################################################################


import numpy as np
import matplotlib.pyplot as plt


##############################################################################


def Plot_Pdelta(F, folder, filename, d_0):
    
    fig, ax = plt.subplots()
    
    ax.set_xlabel('Horizontal displacement [cm]')
    ax.set_ylabel('Force [kN]')
    
    ax.set_title('P-delta curves')
    plt.grid(linestyle='--', linewidth=2)
    
    fname = 'imperfection_0.0' + '/'
    fname += 'Linear' + '/'
    fname += filename
    Data = np.loadtxt(fname, delimiter = ' ')
    Data = np.insert(Data, [0], [[0, 0, 0]], axis= 0)
    
    P = Data[:, 0]
    P *= -F
    P /= 10**3
        
    disp_X = Data[:, 1]
    #disp_X *= -1
    disp_X *= 100
    
    lnstyle = 'dotted'
    
    ax.plot(disp_X, P, label = 'Linear' + ' ' + 'analysis' + ' with d_0 = 0.0', 
            linestyle = lnstyle)
    
    for Analysis in ['Linear', 'Non-Linear']:
        
        fname = folder + '/'
        fname += Analysis + '/'
        fname += filename
        Data = np.loadtxt(fname, delimiter = ' ')
        Data = np.insert(Data, [0], [[0, 0, 0]], axis= 0)
        
        P = Data[:, 0]
        P *= -F
        P /= 10**3
            
        disp_X = Data[:, 1] + d_0
        #disp_X *= -1
        disp_X *= 100
        
        lnstyle = 'solid'
        
        ax.plot(disp_X, P, label = Analysis + ' ' + 'analysis' + ' with d_0 = ' + str(d_0), 
                linestyle = lnstyle)
        
    ax.legend(loc = 'lower right')
    #ax.legend(loc = 'upper left')
    
    #ax.set_xlim(left = 0, right = 1.1)
    #ax.set_ylim(bottom = 0, top = 50)
    
    filename = filename[:-4]
    filename += '.png'
    plt.savefig(filename)
    
    filename = filename[:-4]
    filename += '.pdf'
    plt.savefig(filename)
    
    plt.show()
    

##############################################################################

text1 = '''
Please enter the value of the initial imperfection in cm (if there is none, enter 0):
'''

d_0 = int(input(text1))
d_0 /= 100
folder = 'imperfection_' + str(d_0)


text2 = '''
Please enter the name of the file containing the results of the analysis (don't write the ".txt"):
'''

filename = input(text2)
temp = filename
filename += '.txt'

text3 = '''
Please enter the value of the load that was applied when analizing in kN (P):
'''

P = input(text3)
P = float(P)
P *= 10**3

Plot_Pdelta(P, folder, filename, d_0)


##############################################################################


