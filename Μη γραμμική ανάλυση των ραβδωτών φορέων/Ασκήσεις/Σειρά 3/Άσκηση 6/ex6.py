"""
This program analyzes the structure of excercise 6 taking into account:
    - the initial geometry without any imperfections
    - the material nonlinearities
    and produces the text files with:
    - the vertical displacement of the node where the moment is applied

@author: John Kouretas
"""

import numpy as np
import openseespy.opensees as ops
import opsvis as opsv


##############################################################################


# 1) Define the units:

kg = 1
m  = 1
s  = 1

N = kg * m / s**2
Pa = N / m**2

cm = m / 100
mm = cm / 10

gr  = kg / 10**3
Mg  = kg * 10**3

kN = N * 10**3
MN = kN * 10**3
GN = MN * 10**3

kPa = Pa * 10**3
MPa = kPa * 10**3
GPa = MPa * 10**3

inches = 2.54 * cm
kips   = 4.4482216 * kN
ksi    = 6.89476 * MPa


##############################################################################


# 2) Define auxiliary functions:


    # A) Setup geometry:

def geom(L, h, Nelems):
    
        # Define the coordinates:
    
    x = np.linspace(0, 1, num = Nelems, endpoint=False).reshape((Nelems, 1)) * L
    y = np.linspace(0, 1, num = Nelems, endpoint=False).reshape((Nelems, 1)) * h
    
    x = np.append(x, x + L, axis = 0)
    y = np.append(y, h - y, axis = 0)
    
    coords = np.append(x, y, axis = 1)
    coords = np.append(coords, np.array([[2 * L, 0.00]]), axis = 0)
    
    Nnodes = len(coords)
    #print(coords)
    #print(Nnodes)
    
    ##########################################################################
    
        # Define the nodes:
    
    [ops.node(i + 1, *coords[i, :])    for i in range(Nnodes)]
    
    ##########################################################################
    
    
        # Define boundary conditions:
    
    fix = [[[[0, 0, 0], [0, 1, 0]], [[1, 0, 0], [1, 1, 0]]],
           [[[0, 0, 1], [0, 1, 1]], [[1, 0, 1], [1, 1, 1]]]]
    
            # Fix the first node:
    
    #print(*fix[-1][-1][-1])
    ops.fix(1,      *fix[-1][-1][-1])
    ops.fix(Nnodes, *fix[-1][-1][-1])
    
    ##########################################################################
    
    
    return coords


    # B) Setup elements:

def elements(coords, L, A, E, I, analysis_type = 'Linear'):
    
    
        # Define the number of elements:
    
    Nnodes = len(coords)
    Nelems = Nnodes - 1
    
    
    ##########################################################################
    
    
        # Define the elements:
    
    eleType = 'elasticBeamColumn'
    
    ops.geomTransf('Linear', 1)
    ops.geomTransf('PDelta', 11)
    ops.geomTransf('Corotational', 101)
    
    if analysis_type == 'Linear':
        gTransfTag = 1
        
    elif analysis_type == 'PDelta':
        gTransfTag = 11
        
    else:
        gTransfTag = 101
        
    
    #print(eleType,A,Nelems, gTransfTag)
        # Define the elements:
    
    [ops.element(eleType,
                     i + 1, *[i + 1, i + 2],
                     A, E, I, gTransfTag
                     )   for i in range(Nelems)]
    
    
    ##########################################################################


    # C) Setup recorders:

def recs(coords, analysis_type = 'Linear'):
    
    
    Nnodes = len(coords)
    analysis_Dir = analysis_type + '/'
    
        # Last node's recorder:
    #print(analysis_Dir + 'plotCRV_displ.txt')
    ops.recorder('Node', '-file', analysis_Dir + 'plotCRV_displ.txt',
                     '-time',
                     '-node', *[int((Nnodes + 1) / 2)],
                     '-dof', *[1, 2, 3],
                     'disp')
    
    
##############################################################################


    # D) Setup analysis parameters:

        # P - delta analysis:

def PDelta_analysis(coords, P, u_max, steps, tol, max_iter):
    
    
    ##########################################################################
    
    Nnodes = len(coords)
    
    ##########################################################################
    
    tsTag = 1
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 1
    ops.pattern('Plain', pattTag, tsTag)
    
    ops.load(int((Nnodes + 1) / 2), *P)
    
    ##########################################################################
    
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance', tol, max_iter)
    ops.algorithm('Newton')
    
    temp = u_max / steps
    ops.integrator('DisplacementControl', int((Nnodes + 1) / 2), 2, temp)
    
    ops.analysis('Static')
    
    ops.analyze(steps)
    
    ##########################################################################
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(P, L, h, A, E, I, Nelems, u_max, steps, tol, max_iter,
                                           analysis_type = 'Linear'):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    coords = geom(L, h, Nelems)
    
    opsv.plot_model('nodes', 'elements')
    elements(coords, L, A, E, I, analysis_type)
    
    opsv.plot_model('nodes', 'elements')
    
    recs(coords, analysis_type)
    PDelta_analysis(coords, P, u_max, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 65.715 * cm
L /= 2
h = 0.98 *cm

A = 0.408262 * cm**2

I = 0.0132651 * cm**4

E = 19971.4 * kN / cm**2


    # Define load:

P = np.array([0 * kN, -1 * N, 0 * kN * m])


tol = 1.0e-8
max_iter = 100
steps = 10**3


# 5) :

text = '''
How many elements do you require? (For 7, press enter)
'''

Nelems = input(text)

if len(Nelems) == 0:    Nelems = 7

Nelems = int(Nelems)


text = '''
Please select the type of the analysis:
    For Linear,         press 0.
    For P-Delta,        press 1. (Not recomended) 
    For Corotational,   press 2.
'''

choice = input(text)
choice = int(choice)

if choice == 0:
    analysis_type = 'Linear'
elif choice == 2:
    analysis_type = 'Corotational'
else:
    analysis_type = 'PDelta'

u_max = h
u_max *= 220 / 100
u_max *= P[1] / np.abs(P[1])

##############################################################################


solver(P, L, h, A, E, I, Nelems, u_max, steps, tol, max_iter, analysis_type)


##############################################################################


