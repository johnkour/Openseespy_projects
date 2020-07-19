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
    lst    : LIST OF FLOATS
        Initial list with the values in the Imperial System.
    scalar : FLOAT
        The value used to convert from Imperial to International System.

    Returns
    -------
    result : LIST OF FLOATS
        The list with the values converted to S.I..
        
    '''
    
    result = [t * scalar for t in lst]
    
    return result
    

(x, y, A) = ([t * ft  for t in x], [t * ft  for t in y],
                     [t * inches**2  for t in A])

F = [[t * kip  for t in f]  for f in F]

##############################################################################
#####                     Main Analysis' functions                       #####
##############################################################################

def Model_Build(x, y, A, E):
    '''
    Description
    -----------
    This function is used to determine the basic parameters of the structural
    problem at hand.
    
    Parameters
    ----------
    x : LIST OF FLOATS
        The list of the coordinates of the nodes along the x-axis.
    y : LIST OF FLOATS
        The list of the coordinates of the nodes along the y-axis.
    A : LIST OF FLOATS
        The list with the materials used for the different elements.
    E : FLOAT
        The modulus of elesticity of the elements.
    Returns
    -------
    None.

    '''
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
    
    # Plot model.
    opsplt.plot_model()
    
def Rec_Setup(analysis):
    '''
    Description
    -----------
    This function is used to set up the recorders. It stores the output of the 
    recorders to a folder, whose name is the value of the variable: analysis.

    Parameters
    ----------
    analysis : STRING
        The name of the analysis, currently performed.

    Returns
    -------
    None.

    '''
    
    analysis += '/'
    
    # Set recorders.
    ops.recorder('Node', '-file', analysis + 'NodeDisp.txt', '-time', '-node', *[4], '-dof', *[1, 2, 3], 'disp')
    ops.recorder('Node', '-file', analysis + 'ReactDisp.txt', '-time', '-node', *[4], '-dof', *[1, 2, 3], 'reaction')
    ops.recorder('Element', '-file', analysis + 'ElementsForces.txt', '-time', '-ele', *[1, 2, 3], 'forces')
    

def Analysis_setup(analysis, F, N = 1):
    '''
    Description
    -----------
    This functions is used to setup and then run the analysis.

    Parameters
    ----------
    analysis : STRING
        The name of the analysis, currently performed.
    F        : LIST OF LISTS OF FLOATS
        The list containig a list of loads along the x and y axises and around
        the z axis for every node.
    N        : INTEGER
        The number of the analysises to run. Default value: 1

    Returns
    -------
    None.

    '''
    
    # Auxiliary variable.
    m = len(F)
    
    # Create timeseries.
    ops.timeSeries('Linear', 1)
    
    # Create load pattern.
    ops.pattern('Plain', 1 , 1)

    # Define loads.
    [ops.load(i + 1, *F[i])  for i in range(m)]

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
    
    ok = ops.analyze(N)
    
    if ok == 0:
        status = 'Analysis complete, everything went smoothly.'
    else:
        status = 'ERROR: ANALYSIS FAILED TO CONVERGE' + '\n' + 'Tip: Change algorithm'
    
    print(analysis + '\n' + status)
    
    # Close recorders and scrap model.
    ops.wipe()
    

##############################################################################
#####                            Main Analysis                           #####
##############################################################################

# Step 1: Initilize model parameters.
Model_Build(x, y, A, E)

# Step 2: Name the type of the analysis to be performed.
analysis = 'Static'

# Step 3: Set up the recorders.
Rec_Setup(analysis)

# Step 4: Perform the analysis.
N = 10    # Number of analysises to be performed.
Analysis_setup(analysis, F, N)
