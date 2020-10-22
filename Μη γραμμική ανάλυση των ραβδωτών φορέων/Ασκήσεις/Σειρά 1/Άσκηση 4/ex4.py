"""
This program analyzes the truss of excercise 1 linearly or non-linearly and
produces the text files with:
    - the displacements of the node where the force is applied
    - the internal force of the truss element in the local system
    - the internal force of the truss element in the global system

This time we use the Displacement control instead of the Load control method.

@author: John Kouretas
"""

import numpy as np
import sys
import openseespy.opensees as ops
import openseespy.postprocessing.Get_Rendering as opsplt


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


    # A) Subfunction to define each subsquare of the truss:

def connections(member_type, A, matTag, i, eleTag):
    
    
    ops.element(member_type,     eleTag, *[    i, i + 1], A, matTag)
    ops.element(member_type, eleTag + 1, *[    i, i + 2], A, matTag)
    ops.element(member_type, eleTag + 2, *[    i, i + 3], A, matTag)
    ops.element(member_type, eleTag + 3, *[i + 1, i + 3], A, matTag)
    
    eleTag += 4
    
    return eleTag
    
    
    # B) Setup geometry:

def geom(coords):
    
        # Define the nodes:
    
    [ops.node(i, *coords[i, :])    for i in range(len(coords))]
    
    
    ##########################################################################
    
    
        # Define boundary conditions:
    
    fix = [[1, 1], [1, 0]]
    
            # Pin the first node:
    
    ops.fix(1, *fix[0])
    
            # Roll the second node:
    
    ops.fix(2, *fix[0])
    
    
    ##########################################################################


    # C) Setup elements:

def elements(coords, A, E, analysis_type):
    
    
    matTag = 10
    
    ops.uniaxialMaterial('Elastic', matTag, E)
    
    
    ##########################################################################
    
    
    eleTag = 1
    
    member_type = 'Truss'
    if analysis_type != 'Linear':    member_type = 'corotTruss'
    
    temp = len(coords) - 3
    for i in range(0, temp, 2):
        eleTag = connections(member_type, A, matTag, i, eleTag)
    
    temp += 1
    ops.element(member_type, eleTag, *[temp, temp + 1], A, matTag)
    
    
    ##########################################################################


    # D) Setup recorders:

def recs(coords, analysis_type):
    
    
    temp = len(coords)
    temp -= 1
    
    ops.recorder(   'Node', '-file', analysis_type + '/' +    'plotCRV.txt', '-time', '-node', *[temp],   '-dof', *[2], 'disp')
    
    
##############################################################################


    # E) Setup analysis parameters:

def analysis(coords, P, du, steps, tol, max_iter):
    
    
    tsTag = 1
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 1
    ops.pattern('Plain', pattTag, tsTag)
    
    temp = len(coords) - 1
    ops.load(temp, *P)
    
    
    ##########################################################################
    
    
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance', tol, max_iter)
    ops.algorithm('Newton')
    
    ops.integrator('DisplacementControl', temp, 2, du, max_iter)
    
    ops.analysis('Static')
    
    ops.analyze(steps)
    
    
    ##########################################################################
    
    
    opsplt.plot_model()
    
    
    ##########################################################################
    
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(coords, P, du, steps, tol, max_iter, analysis_type):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    
    geom(coords)
    elements(coords, A, E, analysis_type)
    recs(coords, analysis_type)
    analysis(coords, P, du, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 10 * inches
z = 0.5 * inches

P = np.array([0 * kN, -20 * kips])

A = 0.1 * inches**2
E = 29 * 10**3 * ksi

tol = 1 / 10**6
max_iter = 20
steps = 500
du = -1 * mm


    # Get the type of the analysis from user:

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
    steps = 200
else:
    print('ERROR: THE TYPE OF THE ANALYSIS WAS NOT PROPERLY DEFINED')
    sys.exit()


    # Calculate auxiliary variables:

x = np.linspace(0, 20, num= 21)
x *= z

y = np.linspace(0, 1, num= 2)
y *= z

coords = [[x[i], y[j]]    for i in range(len(x))    for j in range(len(y))]
coords = np.array(coords)


##############################################################################


solver(coords, P, du, steps, tol, max_iter, analysis_type)


##############################################################################


