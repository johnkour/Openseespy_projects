"""
This program analyzes the truss of excercise 1 linearly or non-linearly and
produces the text files with:
    - the displacements of the node where the force is applied
    - the internal force of the truss element in the local system
    - the internal force of the truss element in the global system

Moreover, it provides the user with the choice of calculating the P - delta
curve using the theoretical values from the equation:
    P = (E * A) * (h^2 * u - 1.5 * h * u^2 + 0.5 * u^3) / L^3

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

def geom(coords):
    
        # Define the nodes:
    
    ops.node(1, *coords[0, :])
    ops.node(2, *coords[1, :])
    
    
    ##########################################################################
    
    
        # Define boundary conditions:
    
    fix = [[1, 1], [1, 0]]
    
            # Pin the first node:
    
    ops.fix(1, *fix[0])
    
            # Roll the second node:
    
    ops.fix(2, *fix[1])
    
    
    ##########################################################################


    # B) Setup elements:

def elements(A, E, analysis_type):
    
    
    matTag = 10
    
    ops.uniaxialMaterial('Elastic', matTag, E)
    
    
    ##########################################################################
    
    
    eleTag = 1
    
    member_type = 'Truss'
    if analysis_type != 'Linear':    member_type = 'corotTruss'
    
    ops.element(member_type, eleTag, *[1, 2], A, matTag)
    
    
    ##########################################################################


    # C) Setup recorders:

def recs(analysis_type):
    
    
    ops.recorder(   'Node', '-file', analysis_type + '/' +    'plotCRV.txt', '-time', '-node', *[2],   '-dof', *[2], 'disp')
    ops.recorder('Element', '-file', analysis_type + '/' + 'eleGlobal.txt', '-time',   '-ele', *[1], 'forces')
    ops.recorder('Element', '-file', analysis_type + '/' +  'eleLocal.txt', '-time',   '-ele', *[1], 'basicForces')
    
    
##############################################################################


    # D) Setup analysis parameters:

def analysis(P, steps, tol, max_iter):
    
    
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
    
    
    opsv.plot_model()
    
    
    ##########################################################################
    
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(coords, P, A, E, steps, tol, max_iter, analysis_type):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    
    geom(coords)
    elements(A, E, analysis_type)
    recs(analysis_type)
    analysis(P, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 3 * m
h = 4 * cm

P = np.array([0 * kN, -100 * kN])

A = 1 * m**2
E = 50 * GPa

tol = 1 / 10**6
max_iter = 10
steps = 100


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
    tol = 1 / 10**4
    steps = 200
else:
    print('ERROR: THE TYPE OF THE ANALYSIS WAS NOT PROPERLY DEFINED')
    sys.exit()


    # Calculate auxiliary variables:

temp = np.sqrt(L**2 - h **2)

coords = list()
coords.append([   0, 0])
coords.append([temp, h])

coords = np.array(coords)


##############################################################################


solver(coords, P, A, E, steps, tol, max_iter, analysis_type)


##############################################################################

# 5) Calculate and save the P - delta curve using the theoritical value of P = f(u):

text = '''
Do you want to calculate the P - delta curve using the theoretical solution?
    If  no,     press 0.
    If yes,     press 1.
'''

choice = input(text)
choice = int(choice)

if choice == 1:
    
    steps = 100
    u_max = 3.5 * cm
    
    u = np.arange(steps, dtype = np.float64).reshape([steps, 1])
    u /= steps
    u *= u_max
    
    F = h**2 * u
    F -= 1.5 * h * u**2
    F += 0.5 * u**3
    F *= E * A / L**3
    F /= P[-1]
    
        # Save the curve in *.txt file:
    
    Data = np.concatenate((F, u), axis = 1)
    Data *= -1
    
    path = 'Theoretical/'
    fname = 'plotCRV'
    fname += '.txt'
    
    np.savetxt(path + fname, Data, fmt = '%.8f', delimiter = ' ')

elif choice != 0:
    
    print('ERROR: THE VALUE YOU WROTE WAS OUT OF BOUNDS')
    sys.exit()
