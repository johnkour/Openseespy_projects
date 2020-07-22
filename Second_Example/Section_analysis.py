# Section analysis with Openseespy:

import openseespy.opensees as ops
import openseespy.postprocessing.Get_Rendering as opsplt
import numpy as np

##############################################################################
#####                       Define units in SI                           #####
##############################################################################

# Basic units:
m = 1
kg = 1
s = 1

# Basic structural units:
N = kg * m / s**2
Pa = N / m

kN = 10**3 * N
kPa = 10**3 * Pa

MN = 10**3 * kN
MPa = 10**3 * kPa

GN = 10**3 * MN
GPa = 10**3 * MPa

cm = m / 10**2
mm = cm / 10

##############################################################################
#####                         Input Variables                            #####
##############################################################################

# Model coordinates:
l = 3 * m               # Length of the column.
n = int(l / (1 * m))    # Number of subcolumns - finite elements.

coords = [[0.0, float(i), 0.0]  for i in range(n + 1)]
coords = np.array(coords)
coords *= m

#print(*coords)

# Section geometry:
sec_dims = [[400 * mm, 600 * mm]]
sec_dims = np.array(sec_dims)

#print(sec_dims)

# Rebar coordinates:
rebar_coords_Y = [-225, -150, 150, 225]
rebar_coords_Z = [-150, 0, 150]

rebar_coords_YZ = [[float(rebar_coords_Y[i]), float(rebar_coords_Z[j])]
                                for i in range(len(rebar_coords_Y))
                                for j in range(len(rebar_coords_Z))]

rebar_coords_YZ = np.array(rebar_coords_YZ)
rebar_coords_YZ *= mm

#print(rebar_coords_YZ)

# Concrete properties:
fc = 30.0 * MPa

# Rebar properties:
Es = 200.0 * GPa
fy = 250.0 * MPa
D = 30.0 * mm


##############################################################################
#####                     Main Analysis' functions                       #####
##############################################################################


def get_Sections(fc, Es, fy, sec_dims, D, rebar_coords_YZ):
    '''
    Summary
    -------
    
    
    Parameters
    ----------
    fc : TYPE
        DESCRIPTION.
    Es : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.
    sec_dims : TYPE
        DESCRIPTION.
    D : TYPE
        DESCRIPTION.
    rebar_coords_YZ : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # Define standard materials:
    E_flex = 1.0
    E_rigid = 100.0 * 10**12
    
    ops.uniaxialMaterial('Elastic', 1, Es)
    ops.uniaxialMaterial('Elastic', 10, E_flex)
    ops.uniaxialMaterial('Elastic', 20, E_rigid)
    
    # Rebar material:
    
    b = 0.02
    params = [18.0, 0.925, 0.15]   # Recommended values for the Menegotto-Pinto model.
    #   uniaxialMaterial('Steel02', matTag, Fy, E0, b, *params)
    ops.uniaxialMaterial('Steel02',    101, fy, Es, b, *params)
    
    # Concrete material:
    
    fc *= -1                # Value from concrete lab, NTUA.
    fcu = 0.2 * fc          # Value from concrete lab, NTUA.
    eps_c0 = -2 / 10**3     # Value from concrete lab, NTUA.
    eps_u = -3.5 / 10**3    # Value from concrete lab, NTUA.
    
    #   uniaxialMaterial('Concrete01', matTag, fpc,  epsc0, fpcu,  epsU)
    ops.uniaxialMaterial('Concrete01',    201,  fc, eps_c0,  fcu, eps_u)
    
    # Geometry preprocessing:
    vertices = np.zeros((4,2))
    vertices[:2, :] = -sec_dims / 2
    vertices[2:, :] = sec_dims / 2
    vertices[0, 0] *= -1
    vertices[2, 0] *= -1
    
    N_fiber_Y = 50
    N_fiber_Z = 1
    
    A_bar = np.pi * D**2 /4
    
    #print(vertices)
    # Define sections:
    #   For the rebar/concrete part we use the fiber/patch command
        # Define test sections:
    
    #   section('Fiber', secTag)
    ops.section('Fiber',      1)
    
    #   patch('quad', matTag, numSubdivIJ, numSubdivJK,          *crdsI,          *crdsJ,          *crdsK,          *crdsL)
    ops.patch('quad',      1,   N_fiber_Z,   N_fiber_Y, *vertices[0, :], *vertices[1, :], *vertices[2, :], *vertices[3, :])
    
    n = rebar_coords_YZ.shape[0]
    #    fiber(                 yloc,                  zloc,     A, matTag)
    #ops.fiber(rebar_coords_YZ[i, 0], rebar_coords_YZ[i, 1], A_bar,      1)
    [ops.fiber(rebar_coords_YZ[i, 0], rebar_coords_YZ[i, 1], A_bar,      1)
                             for i in range(n)]
    
        # Define actual sections:
    
    #   section('Fiber', secTag)
    ops.section('Fiber',     10)
    
    #   patch('quad', matTag, numSubdivIJ, numSubdivJK,          *crdsI,          *crdsJ,          *crdsK,          *crdsL)
    ops.patch('quad',    201,   N_fiber_Z,   N_fiber_Y, *vertices[0, :], *vertices[1, :], *vertices[2, :], *vertices[3, :])
    
    n = rebar_coords_YZ.shape[0]
    #    fiber(                 yloc,                  zloc,     A, matTag)
    #ops.fiber(rebar_coords_YZ[i, 0], rebar_coords_YZ[i, 1], A_bar,    101)
    [ops.fiber(rebar_coords_YZ[i, 0], rebar_coords_YZ[i, 1], A_bar,    101)
                             for i in range(n)]
    
    # Geometric tranformation:
    #   geomTransf(transfType, transfTag, *transfArgs)
    ops.geomTransf(  'Linear',         1)             # test
    ops.geomTransf(  'PDelta',        10)             # actual
    
    # Integrator:
    N = 4
    #   beamIntegration('Lobatto', tag, secTag, N)
    ops.beamIntegration('Lobatto',   1,      1, N)    # test
    ops.beamIntegration('Lobatto',  10,     10, N)    # actual


def get_Model(coords):
    '''
    Summary
    -------
    
    
    Parameters
    ----------
    coords : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # Define nodes:
    N = coords.shape[0]
    [ops.node(i + 1, *coords[i, :])  for i in range(N)]
    
    # Fix base point:
    fixxity = [1, 1, 1]         # The base is fully fixxed.
    ops.fix(1, *fixxity)
    
    # Define elements:
    transfTag = 1       # testing.
    #transfTag = 10      # actual.
    integrationTag = 1  # testing.
    #integrationTag = 10 # actual.
    #    element('forceBeamColumn', eleTag,   *eleNodes, transfTag, integrationTag, '-iter', maxIter=10, tol=1e-12)
    #ops.element('forceBeamColumn',      i, *[i, i + 1], transfTag, integrationTag, '-iter',         30,     1e-12)
    [ops.element('forceBeamColumn',      i, *[i, i + 1], transfTag, integrationTag, '-iter',         30,     1e-12)
                             for i in range(1, N)]
    
    print('Model built.')


##############################################################################
#####                            Main Analysis                           #####
##############################################################################

# Initialization:
ops.wipe()

# Modelbuild:
ops.model('Basic', '-ndm', 2, '-ndf', 3)
get_Sections(fc, Es, fy, sec_dims, D, rebar_coords_YZ)
get_Model(coords)
opsplt.plot_model()

# Analysis:


# Close:
ops.wipe()
