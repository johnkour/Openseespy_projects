"""
This program analyzes the structure of excercise 5 taking into account:
    - the initial geometry without any imperfections
    - the geometric nonlinearities (by assuming a small imperfection)
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

def geom(L, Nelems):
    
        # Define the coordinates:
    
    coords = np.linspace(0, L, Nelems + 1).reshape((Nelems + 1, 1))
    coords = np.append(coords, np.zeros((Nelems + 1, 1)), axis = 1)
    
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
    
    ops.fix(1, *fix[1][1][1])
    
    
    ##########################################################################
    
    
    return coords


    # B) Setup elements:

def elements(coords, L, A, E, Iz, analysis_type = 'Linear'):
    
    
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
        
    
    #print(eleType,A,matTag,Nelems)
        # Define the elements:
    
    [ops.element(eleType,
                     i + 1, *[i + 1, i + 2],
                     A, E, Iz, gTransfTag
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
                     '-node', *[int(Nnodes - 1)],
                     '-dof', *[1, 2, 3],
                     'disp')
    
    
##############################################################################


    # D) Setup analysis parameters:

        # P - delta analysis:

def PDelta_analysis(coords, P, steps, tol, max_iter):
    
    
    ##########################################################################
    
    Nnodes = len(coords)
    
    ##########################################################################
    
    tsTag = 1
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 1
    ops.pattern('Plain', pattTag, tsTag)
    
    ops.load(Nnodes, *P)
    
    ##########################################################################
    
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance', tol, max_iter)
    ops.algorithm('Newton')
    
    temp = 1 / steps
    ops.integrator('LoadControl', temp)
    
    ops.analysis('Static')
    
    ops.analyze(steps)
    
    ##########################################################################
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(P, L, A, E, Iz, Nelems, steps, tol, max_iter, analysis_type = 'Linear'):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    coords = geom(L, Nelems)
    elements(coords, L, A, E, Iz, analysis_type)
    
    opsv.plot_model('nodes', 'elements')
    
    recs(coords, analysis_type)
    PDelta_analysis(coords, P, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 10 * inches

A = 4 * inches**2

Iz = 1.33 * inches**4

E = 100 * ksi


    # Define load:

Mc = np.pi * E * Iz
Mc /= L

P = np.array([0 * kips, 0 * kips, Mc])


tol = 1.0e-6
max_iter = 100
steps = 10000


# 5) :

text = '''
How many elements do you require? (For 17, press enter)
'''

Nelems = input(text)

if len(Nelems) == 0:    Nelems = 17

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

##############################################################################


solver(P, L, A, E, Iz, Nelems, steps, tol, max_iter, analysis_type)


##############################################################################


