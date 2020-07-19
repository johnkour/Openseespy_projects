# Basic truss example in Openseespy:

import openseespy.opensees as ops
import openseespy.postprocessing.Get_Rendering as opsplt

##############################################################################
#####                       Define units in SI                           #####
##############################################################################

# Basic units:
m = 1
kg = 1
s = 1

N = kg * m / s**2
Pa = N / m

inches = 0.0254 * m
ft = 12 * inches
kip = 4.45 * 10**3 * N
ksi = 6.89 * 10**6 * Pa

##############################################################################
#####                         Input Variables                            #####
##############################################################################

x = [0.0, 12.0, 14.0, 6.0]
# x = list(map(lambda a: a * ft, x))

y = [0.0, 0.0, 0.0, 8.0]
# y = list(map(lambda a: a * ft, y))

A = [10.0, 5.0]

E = 3 * 10**3 * ksi

F = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [100.0, -50.0, 0.0]]

def aux_function(lst, scalar):
    '''
    
    Description
    -----------
    This function is used multiply each value in a list with a scalar number.
    It will be very usefull when converting the variables from the Imperial
    system to S.I.

    Parameters
    ----------
    lst    : TYPE
        DESCRIPTION.
    scalar : TYPE
        DESCRIPTION.

    Returns
    -------
    result : TYPE
        DESCRIPTION.
        
    '''
    
    result = [t * scalar for t in lst]
    
    return result
    

(x, y, A) = ([t * ft  for t in x], [t * ft  for t in y],
                     [t * inches**2  for t in A])

F = [[t * kip  for t in f]  for f in F]

##############################################################################
#####                          Main Analysis                             #####
##############################################################################

# Delete existing model.
ops.wipe()

# Define the model.
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Define materials.
ops.uniaxialMaterial('Elastic', 1, E)

# Define the nodes.
m = len(x)

[ops.node(i + 1, *[x[i], y[i]])  for i in range(m)]

# Fix the nodes.
fixxity = [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]]

[ops.fix(i + 1, *fixxity[3])  if i + 1 != 4  else ops.fix(i + 1, *fixxity[0])
             for i in range(m)]

# Define elements.
conn = [[1, 4], [2, 4], [3, 4]]
[ops.element('Truss', i + 1, *conn[i], A[1], 1)  if i != 0
             else  ops.element('Truss', i + 1, *conn[i], A[0], 1)
             for i in range(len(conn))]

# Create timeseries.
ops.timeSeries('Linear', 1)

# Create load pattern.
ops.pattern('Plain', 1 , 1)

# Define loads.
[ops.load(i + 1, *F[i])  for i in range(m)]

# Set recorders.
ops.recorder('Node', '-file', 'NodeDisp.txt', '-time', '-node', *[4], '-dof', *[1, 2, 3], 'disp')
ops.recorder('Node', '-file', 'ReactDisp.txt', '-time', '-node', *[4], '-dof', *[1, 2, 3], 'reaction')
ops.recorder('Element', '-file', 'ElementsForces.txt', '-time', '-ele', *[1, 2, 3], 'forces')

# Define system.
ops.system('BandSPD')

# Define numberer.
ops.numberer('RCM')

# Define constraint handler
ops.constraints('Plain')

# Define integrator.
ops.integrator('LoadControl', 1.0)

# Define algorithm
ops.algorithm('Linear')

# Create analysis object
ops.analysis('Static')

# Execute the analysis
ops.initialize()    # Set recorders to start recording at 0 time.

ok = ops.analyze(10)

if ok == 0:
    status = 'Analysis complete, everything went smoothly.'
else:
    status = 'ERROR: ANALYSIS FAILED TO CONVERGE' + '\n' + 'Tip: Change algorithm'

print(status)

# Plot model.
opsplt.plot_model()

# Close recorders and scrap model.
ops.wipe()
