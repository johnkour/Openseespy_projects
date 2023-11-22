"""
This program analyzes the structure of excercise 3 taking into account:
    - the initial geometry without any imperfections
    - the geometric nonlinearities (by assuming a small imperfection)
    and produces the text files with:
    - the horizontal displacement of the node where the force is applied

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


##############################################################################


# 2) Define auxiliary functions:


    # A) Setup geometry:

def geom(L, d_0 = 0.00):
    
        # Define the coordinates:
    
    coords = np.array([[0         , L[1]],
                       [L[0] + d_0, L[1]],
                       [L[0]      ,    0]])
    
    Nnodes = len(coords)
    #print(coords)
    
    ##########################################################################
    
        # Define the nodes:
    
    [ops.node(i + 1, *coords[i, :])    for i in range(Nnodes)]
    
    
    ##########################################################################
    
    
        # Define boundary conditions:
    
    fix = [[1, 1], [1, 0], [0, 1]]
    
            # Fix the first node:
    
    ops.fix(1, *fix[0])
    ops.fix(3, *fix[0])
    
    
    ##########################################################################
    
    
    return coords


    # B) Setup elements:

def elements(coords, L, A, E, analysis_type = 'Linear'):
    
    
        # Define the number of elements:
    
    Nnodes = len(coords)
    Nelems = Nnodes - 1
    
    
    ##########################################################################
    
    
        # Define elastic material:
    
    matTag = 10
    
    ops.uniaxialMaterial('Elastic', matTag, E)
    
    
    ##########################################################################
    
    
    eleType = 'Truss'
    if analysis_type != 'Linear':   eleType = 'corotTruss'
    
    #print(eleType,A,matTag,Nelems)
        # Define the elements:
    
    [ops.element(eleType,
                     i + 1, *[i + 1, i + 2],
                     A[i], matTag,
                     )   for i in range(Nelems)]
    
    
    ##########################################################################


    # C) Setup recorders:

def recs(coords, d_0 = 0.00, analysis_type = 'Linear'):
    
    
    Nnodes = len(coords)
    analysis_Dir = 'imperfection_' + str(d_0)[:4] + '/'
    analysis_Dir += analysis_type + '/'
    
        # Last node's recorder:
    #print(analysis_Dir + 'plotCRV_displ.txt')
    ops.recorder('Node', '-file', analysis_Dir + 'plotCRV_displ.txt',
                     '-time',
                     '-node', *[int(Nnodes - 1)],
                     '-dof', *[1, 2],
                     'disp')
    
    
##############################################################################


    # D) Setup analysis parameters:

        # P - delta analysis:

def PDelta_analysis(coords, P, steps, tol, max_iter):
    
    
    ##########################################################################
    
    #Nnodes = len(coords)
    
    ##########################################################################
    
    tsTag = 1
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 1
    ops.pattern('Plain', pattTag, tsTag)
    
    ops.load(2, *P)
    
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

def solver(P, L, A, E, d_0, steps, tol, max_iter, analysis_type = 'Linear'):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    
    coords = geom(L, d_0)
    elements(coords, L, A, E, analysis_type)
    
    opsv.plot_model('nodes', 'elements')
    
    recs(coords, d_0, analysis_type)
    PDelta_analysis(coords, P, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = list()

Lab = 4 * m
L.append(Lab)

Lbc = 4 * m
L.append(Lbc)


A = list()

Aab = 2 * mm**2
A.append(Aab)

Abc = 5 * 10**3 * mm**2
A.append(Abc)

E = 200 * GPa


    # Define load using a stiffness parameter:


P = np.array([0 * kN, -400 * kN])


tol = 1.0e-6
max_iter = 100
steps = 10000


    # Set auxiliary variables:

#d_0 = 5 * cm


# 5) Calculate and save the P - delta curve using the theoritical value of P = f(u):

text = '''
Do you want to set a horizontal imperfection of 5cm at node 2?
    If yes,     press 1.
    If  no,     press 0.
'''

choice = input(text)
choice = int(choice)

if choice == 1:
    d_0 = 5 * cm
elif choice == 0:
    d_0 = 0 * cm
else:
    print('ERROR: YOUR ANSWER IS OUT OF BOUNDS')
    sys.exit()

#print(d_0)

text = '''
Please select the type of the analysis:
    For Linear,     press 0.
    For Non-Linear, press 1.
'''

choice = input(text)
choice = int(choice)

if choice == 0:
    analysis_type = 'Linear'
elif choice == 1:
    analysis_type = 'Non-Linear'
else:
    print('ERROR: THE TYPE OF THE ANALYSIS WAS NOT PROPERLY DEFINED')
    sys.exit()

##############################################################################


solver(P, L, A, E, d_0, steps, tol, max_iter, analysis_type)


##############################################################################


