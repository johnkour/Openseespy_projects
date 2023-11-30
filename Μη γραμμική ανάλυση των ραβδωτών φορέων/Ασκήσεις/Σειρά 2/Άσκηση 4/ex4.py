"""
This program analyzes the cantilevel of excercise 4 taking into account the
geometric nonlinearities and produces the text files with:
    - the displacements of the node where the force is applied
    - the basic forces of the cantilevel elements
    - the internal forces of the cantilevel elements in the local system
    - the internal forces of the cantilevel elements in the global system
    - the basic displacements of the cantilevel elements

Moreover, a force of 50kN is applied at the edge of the cantilevel and the
appropriate text files are once again produced.

Lastly, the analysis is repeated for several different numbers of finite
elements (1, 3, 5, 10, 15, 20, 30, 50).

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


##############################################################################


# 2) Define auxiliary functions:


    # A) Setup geometry:

def geom(L, Nnodes):
    
        # Define the coordinates:
    
    coords = np.linspace(0, L, Nnodes)
    coords = coords.reshape((Nnodes, 1))
    
    temp = np.zeros((Nnodes, 1))
    
    coords = np.append(coords, temp, axis= 1)
    
    
    ##########################################################################
    
        # Define the nodes:
    
    [ops.node(i + 1, *coords[i, :])    for i in range(Nnodes)]
    
    
    ##########################################################################
    
    
        # Define boundary conditions:
    
    fix = [[1, 1, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0]]
    
            # Fix the first node:
    
    ops.fix(1, *fix[0])
    
    
    ##########################################################################


    # B) Setup elements:

def elements(Nnodes, L, A, E, Iz):
    
    
        # Define the number of elements:
    
    Nelems = Nnodes - 1
    
    
    ##########################################################################
    
    
        # Define the geometric nonlinearity:
    
    TransfTag = 1
    
    ops.geomTransf('Corotational', TransfTag)
    
    
    ##########################################################################
    
    
        # Define the elements:
    
    [ops.element('elasticBeamColumn',
                     i + 1, *[i + 1, i + 2],
                     A, E, Iz,
                     TransfTag
                     )   for i in range(Nelems)]
    
    
    ##########################################################################


    # C) Setup recorders:

def recs(Nnodes, Nnodes_dir):
    
    
    Nelems = Nnodes - 1
    
        # Last node's recorder:
    
    ops.recorder('Node', '-file', Nnodes_dir + 'plotCantilevel_displ.txt',
                     '-time',
                     '-node', *[Nnodes],
                     '-dof', *[2],
                     'disp'
                     )
    
    
        # Elements' recorders:
    
    ops.recorder('Element', '-file', Nnodes_dir + 'eleBasicForce.txt',
                     '-time',
                     '-ele', *[i + 1    for i in range(Nelems)],
                     'basicForces'
                     )
    
    ops.recorder('Element', '-file',  Nnodes_dir + 'eleLocalForces.txt',
                     '-time',
                     '-ele', *[i + 1    for i in range(Nelems)],
                     'localForce'
                     )
    
    ops.recorder('Element', '-file', Nnodes_dir + 'eleGlobalForces.txt',
                     '-time',
                     '-ele', *[i + 1    for i in range(Nelems)],
                     'globalForce'
                     )
    
    ops.recorder('Element', '-file', Nnodes_dir + 'eleBasicDispl.txt',
                     '-time',
                     '-ele', *[i + 1    for i in range(Nelems)],
                     'basicDeformation'
                     )
    
    
##############################################################################


    # D) Setup analysis parameters:

        # Point load analysis:

def PointLoad_analysis(Nnodes, P, tol, max_iter):
    
    
    tsTag = 10
    ops.timeSeries('Linear', tsTag)
    
    pattTag = 10
    ops.pattern('Plain', pattTag, tsTag)
    
    ops.load(Nnodes, *P)
    
    ##########################################################################
    
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormUnbalance', tol, max_iter)
    ops.algorithm('Newton')
    
    ops.integrator('LoadControl', 1)
    
    ops.analysis('Static')
    
    ops.analyze(1)
    
    ##########################################################################
    
    ops.wipe()
    
    
    ##########################################################################
    
    
        # P-delta analysis:
    
def PDelta_analysis(Nnodes, P, steps, tol, max_iter):
    
    
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
    
    
    opsv.plot_model('nodes')
    
    
    ##########################################################################
    
    
    ops.wipe()
    
    
##############################################################################


# 3) Define the solver:

def solver(Nnodes, P, F, L, A, E, Iz, steps, tol, max_iter, Nnodes_dir):
    
    
    ops.wipe()
    
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    geom(L, Nnodes)
    elements(Nnodes, L, A, E, Iz)
    
    temp = Nnodes_dir + 'PointLoad/'
    
    recs(Nnodes, temp)
    PointLoad_analysis(Nnodes, F, tol, max_iter)
    
    
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    geom(L, Nnodes)
    elements(Nnodes, L, A, E, Iz)
    
    temp = Nnodes_dir + 'PDelta/'
    
    recs(Nnodes, temp)
    PDelta_analysis(Nnodes, P, steps, tol, max_iter)


##############################################################################



# 4) Define the variables of the problem:

L = 4 * m
A = 42 * cm**2
Iz = 6482 * cm**4
E = 200 * GPa


    # Define load using a stiffness parameter:

k = E * Iz
k /= L**2

P = np.array([0 * kN, -10 * k, 0 * kN * m])

F = np.array([0 * kN, -50 * kN, 0 * kN * m])


tol = 1.0e-3
max_iter = 15
steps = 50


    # Set auxiliary variables:

Nnodes = [2, 4, 11, 16, 21, 31, 51]


##############################################################################


for nnodes in Nnodes:
    
    Nnodes_dir = str(nnodes - 1)
    Nnodes_dir += 'FiniteElements/'
    solver(nnodes, P, F, L, A, E, Iz, steps, tol, max_iter, Nnodes_dir)


##############################################################################


