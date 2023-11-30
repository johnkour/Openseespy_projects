"""
A basic program to plot in the same PDF file P - delta curves from text files 
generated from Openseespy recorders.

@author: John Kouretas
"""


##############################################################################


import numpy as np
import sys
import matplotlib.pyplot as plt


##############################################################################


def Plot_Pdelta(F, folders, filename, Nelements, choise):
    
    fig, ax = plt.subplots()
    
    ax.set_xlabel('Displacement [cm]')
    ax.set_ylabel('Force [kN]')
    
    ax.set_title(choise)
    plt.grid(linestyle='--', linewidth=2)
    
    mrk = ''
    
    for i in range(len(folders)):
        
        #if Nelements[i] < 10: continue
        
        fname = folders[i]
        fname += '/' + filename
        Data = np.loadtxt(fname, delimiter = ' ')
        #print(len(Data))
        #print(Data.shape)
        if Data.shape == (len(Data),):
            Data = Data.reshape((1, len(Data)))
            mrk = 'o'
            lc = 'lower right'
        else:
            Data = np.insert(Data, [0], [[0, 0]], axis= 0)
            lc = 'upper left'
        
        #print(Data)
        P = Data[:, 0]
        P *= -F
        P /= 10**3
            
        disp = Data[:, 1]
        disp *= -1
        disp *= 100
        
        ax.plot(disp, P, label = '{:d} Elements'.format(Nelements[i]),
                                        marker = mrk)
        
    ax.legend(loc = lc)
    
    
    #ax.set_xlim(left = 0, right = 1.1)
    #ax.set_ylim(bottom = 0, top = 50)
    
    plt.savefig(choise + '.png')
    plt.savefig(choise + '.pdf')
    
    plt.show()
    

##############################################################################

text = '''
Do you want to compare P - Delta curves or Single Load Cases?
-- For P - Delta curves, press 0.
-- For Single Load Cases, press 1.
'''

choise = input(text)
choise_2 = ''

if choise == '0':
    choise = 'PDelta'
    choise_2 = 'P - Delta curves'
    
elif choise == '1':
    choise = 'PointLoad'
    choise_2 = 'Single Load Cases'
    
else:
    print('ERROR: THE VALUE YOU WROTE WAS OUT OF BOUNDS')
    sys.exit()

Nelements = np.array([1, 3, 10, 15, 20, 30, 50])

folders = list()

[folders.append("{:d}FiniteElements".format(Nelements[i]) + "/" + choise)
                                         for i in range(len(Nelements))]

filename = 'plotCantilevel_displ'
filename += '.txt'

text2 = '''
Please enter the value of the load that was applied when analizing in kN (P):
'''

P = input(text2)
P = float(P)
P *= 10**3

Plot_Pdelta(P, folders, filename, Nelements, choise_2)


##############################################################################


