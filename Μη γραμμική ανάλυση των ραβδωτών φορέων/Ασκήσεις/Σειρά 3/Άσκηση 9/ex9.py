"""
This program analyzes the structure of excercise 8 taking into account:
    - the initial geometry without any imperfections
    - the material nonlinearities
    and produces the text files with:
    - the vertical displacement of the node where the moment is applied

@author: John Kouretas
"""

import numpy as np
import sys
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

def geom(L, Nelems):
    
        # Define the coordinates:
    
    a = np.linspace(0, 1, num = Nelems + 1, endpoint = True).reshape((Nelems + 1, 1)) * L
    b = np.ones((Nelems, 1)) * L
    
    x = np.append(b * 0,    a, axis = 0)
    y = np.append(a,        b, axis = 0)
    
    coords = np.append(x, y, axis = 1)
    
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
    ops.fix(1,          *fix[0][-1][-1])
    ops.fix(Nnodes,     *fix[0][-1][-1])
    
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
    ops.recorder('Node', '-file', analysis_Dir + "Lee's Frame.txt",
                     '-time',
                     '-node', *[int((Nnodes + 1) / 2 + 0.2 * (Nnodes - 1) / 2)],
                     '-dof', *[1, 2, 3],
                     'disp')
    
    #print(int((Nnodes + 1) / 2 + 0.2 * (Nnodes - 1) / 2))
    
##############################################################################


    # D) Setup analysis parameters:

        # P - delta analysis:

def PDelta_analysis(coords, P, steps, ArcLngthParams, tol, max_iter):
    
    
    ##########################################################################
    
    Nnodes = len(coords)
    
    ds = ArcLngthParams[0]
    alpha = ArcLngthParams[1]
    
    ##########################################################################
    
    tsTag = 1
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 1
    ops.pattern('Plain', pattTag, tsTag)
    
    ops.load(int((Nnodes + 1) / 2 + 0.2 * (Nnodes - 1) / 2), *P)
    
    ##########################################################################
    
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', tol, max_iter)
    #ops.test('NormUnbalance', tol, max_iter)
    ops.algorithm('BFGS')
    
    ops.integrator('ArcLength', ds, alpha)
    
    ops.analysis('Static')
    
    ops.analyze(steps)
    
    ##########################################################################
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(P, L, A, E, I, Nelems, steps, ArcLngthParams, tol, max_iter,
                                           analysis_type = 'Linear'):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    coords = geom(L, Nelems)
    
    elements(coords, L, A, E, I, analysis_type)
    
    opsv.plot_model('nodes', 'elements')

    recs(coords, analysis_type)
    PDelta_analysis(coords, P, steps, ArcLngthParams, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 1.20 * m

A = 6 * cm**2

I = 2 * cm**4

E = 7060.8 * kN / cm**2


    # Define load:

P = np.array([0 * N, -100 * kN, 0 * kN * m])

tol = 1.0e-5
max_iter = 200
steps = 10**6


# 5) :

text = '''
How many elements do you require? (For 10, press enter)
'''

Nelems = input(text)

if len(Nelems) == 0:    Nelems = 10

Nelems = int(Nelems)

if Nelems == 1:
    print('ERROR: THE TYPE OF THE ANALYSIS DOES NOT SUPPORT ONE (1) ELEMENT')
    sys.exit()


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

ArcLngthParams = list()


s = np.abs(P)[1] / 10**2 * 3
ds = s / steps
alpha = 1

ArcLngthParams.append(ds)
ArcLngthParams.append(alpha)

##############################################################################


solver(P, L, A, E, I, Nelems, steps, ArcLngthParams, tol, max_iter,
                                                               analysis_type)


##############################################################################


