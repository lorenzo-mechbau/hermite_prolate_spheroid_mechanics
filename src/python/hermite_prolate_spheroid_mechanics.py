#!/usr/bin/env python
import os
from numpy import pi
from opencmiss.iron import iron

import prolate_spheroid_geometry


# Prolate spheroid geometry parameters:
cutoffAngle = 120.0 * pi / 180.0
focus = 37.5  # mm
endocardiumLambda = 0.38
epicardiumLambda = 0.69

# Fibre angles in radians:
epicardiumFibreAngle = 45.0 * pi / 180.0
endocardiumFibreAngle = -90.0 * pi / 180.0
# Assume constant sheet angle
sheetAngle = 90.0 * pi / 180.0

# Simulation parameters:
cavityPressure = 0.8
numIncrements = 1

# Geometric and hydrostatic pressure interpolation:
interpolations = ['cubic_hermite', 'linear']
# Other interpolations should also work, eg:
#interpolations = ['quadratic', 'linear']
nodalPressure = True
geometricMeshComponent = 1
pressureMeshComponent = 2
hasDerivatives = interpolations[0] == 'cubic_hermite'
# Circumferential, longitudinal, transmural number of elements:
numberGlobalElements = [4, 2, 1]

# Set up a ProlateSpheroid object, which calculates the mesh geometry:
geometry = prolate_spheroid_geometry.ProlateSpheroid(
        focus, endocardiumLambda, epicardiumLambda, cutoffAngle,
        numberGlobalElements, endocardiumFibreAngle, epicardiumFibreAngle, sheetAngle,
        interpolations)

# Constitutive relation setup
# Guccione constitutive relation:
constitutiveRelation = iron.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_EXPONENTIAL
constitutiveParameters = [0.88, 0.0, 18.5, 3.58, 3.26]
initialHydrostaticPressure = 0.0

# User numbers for identifying OpenCMISS objects
coordinateSystemUserNumber = 1
regionUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
(geometricFieldUserNumber,
    fibreFieldUserNumber,
    materialFieldUserNumber,
    dependentFieldUserNumber,
    deformedFieldUserNumber,
    equationsSetFieldUserNumber) = range(1, 7)
equationsSetUserNumber = 1
problemUserNumber = 1

# Get the number of computational nodes and this computational node number
#computationEnvironment = iron.ComputationEnvironment()
#numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
#computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()


# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region within the world region and
# assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.LabelSet("ProlateSpheroid")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Create mesh from the prolate spheroid geometry
mesh = geometry.generateMesh(region)

# Create a decomposition for the mesh
# This breaks the mesh up into multiple decomposition
# domains when running in parallel
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
# Have to enable face calculation for the decomposition
# in order to be able to use pressure boundary conditions
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
# Set the x, y and z components to use the first mesh component:
for component in range(1, 4):
    geometricField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, component, geometricMeshComponent)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
geometricField.CreateFinish()

# Use the prolate spheroid geometry to set the geometric field parameters
#geometry.setGeometry(computationEnvironment,geometricField)
geometry.setGeometry(geometricField)

# Create a fibre field and attach it to the geometric field
# This has three components describing fibre orientations as
# rotations in radians about the z, y' and x'' base vectors.
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U, "Fibre")
fibreField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
# Use linear mesh component for fibre field
for component in range(1, 4):
    fibreField.ComponentMeshComponentSet(
                iron.FieldVariableTypes.U, component, pressureMeshComponent)
fibreField.CreateFinish()

# Use the prolate spheroid geometry to set up the fibre field values
#geometry.setFibres(computationEnvironment,fibreField)
geometry.setFibres(fibreField)

# Create the equations set and equations set field
# This defines the type of equations to solve
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY, 
                             iron.EquationsSetTypes.FINITE_ELASTICITY,
                             constitutiveRelation]
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the material field, used for setting constitutive parameters
# This has an interpolation type of constant by default
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

# Set constant material parameters:
for (component, value) in enumerate(constitutiveParameters, 1):
    materialField.ComponentValuesInitialise(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
            component, value)

# Create the dependent field for storing the solution
# This has one U variable for displacement values and hydrostatic pressure,
# and one DELUDELN variable for the nodal forces
# Sensible defaults are defined by the equations set, but we can
# override some settings here, eg. to change the interpolation type
# for the hydrostatic pressure to node based
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
if nodalPressure:
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4,
            iron.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4,
            iron.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, 4, pressureMeshComponent)
    dependentField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.DELUDELN, 4, pressureMeshComponent)
else:
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.U, 4,
            iron.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(
            iron.FieldVariableTypes.DELUDELN, 4,
            iron.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
equationsSet.DependentCreateFinish()

# Initialise dependent field position values from the undeformed geometry
for component in [1, 2, 3]:
    geometricField.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component,
        dependentField, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component)
    dependentField.ComponentValuesInitialise(
        iron.FieldVariableTypes.DELUDELN,
        iron.FieldParameterSetTypes.VALUES,
        component, 0.0)

# Set the initial hydrostatic pressure
dependentField.ComponentValuesInitialise(
    iron.FieldVariableTypes.U,
    iron.FieldParameterSetTypes.VALUES,
    4, initialHydrostaticPressure)
dependentField.ComponentValuesInitialise(
    iron.FieldVariableTypes.DELUDELN,
    iron.FieldParameterSetTypes.VALUES,
    4, 0.0)

# Create a deformed geometry field, as cmgui doesn't like displaying
# deformed fibres from the dependent field because it isn't a geometric field.
deformedField = iron.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.MeshDecompositionSet(decomposition)
deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
for component in [1, 2, 3]:
    deformedField.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, component,
            geometricMeshComponent)
deformedField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
deformedField.CreateFinish()

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Define the problem
# The problem defines how the equations are solved by
# setting up control loops and solvers
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.NONE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create the problem control loops
# For static finite elasticity, there is just a
# single load increment control loop, and we set
# the number of load increments by setting the number of
# iterations of this control loop
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.MaximumIterationsSet(numIncrements)
controlLoop.OutputTypeSet(iron.ControlLoopOutputTypes.PROGRESS)
problem.ControlLoopCreateFinish()

# Create the problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.OutputTypeSet(iron.SolverOutputTypes.MONITOR)
solver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
solver.NewtonRelativeToleranceSet(1.0e-14)
solver.NewtonAbsoluteToleranceSet(1.0e-14)
solver.NewtonSolutionToleranceSet(1.0e-14)
# Adjust settings for the line search solver
linesearchSolver = iron.Solver()
solver.NewtonLinearSolverGet(linesearchSolver)
linesearchSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

# Create solver equations for the problem and add
# the single equations set to solver equations.
# For coupled problems there may be multiple equations sets
# solved within one solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
# Only a single equations set to add, the finite elasticity equations:
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions on the solver equations
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

def getDomainNodes(geometry, decomposition, component):
#def getDomainNodes(computationEnvironment, geometry, decomposition, component):
    component_name = interpolations[component - 1]
    computationalNodeNumber = iron.ComputationalNodeNumberGet()
#    computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
    nodes = geometry.componentNodes(component_name)
    meshComponent = geometry.meshComponent(component_name)
    return set(node for node in nodes
        if decomposition.NodeDomainGet(node, meshComponent) == computationalNodeNumber)
geometricDomainNodes = getDomainNodes(geometry, decomposition, geometricMeshComponent) 
# geometricDomainNodes = getDomainNodes(computationEnvironment, geometry, decomposition, geometricMeshComponent)

# Fix epicardium nodes at the base:
baseNodes = set(geometry.nodeGroup('base'))
externalNodes = set(geometry.nodeGroup('external'))
fixedNodes = baseNodes.intersection(externalNodes)
for node in fixedNodes.intersection(geometricDomainNodes):
    for component in (1, 2, 3):
        derivative, version = 1, 1
        # AddNode is used, as this fixes the nodal value at an
        # increment (in this case 0) from the current value
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U,
                version, derivative, node, component,
                iron.BoundaryConditionsTypes.FIXED, 0.0)
        if hasDerivatives:
            derivative = iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1
            boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U,
                    version, derivative, node, component,
                    iron.BoundaryConditionsTypes.FIXED, 0.0)

# Constrain degrees of freedom at the apex to collapse faces:
constrainedNodeSets = geometry.constrainedNodes()
for nodes in constrainedNodeSets:
    in_domain = [n in geometricDomainNodes for n in nodes]
    if all(in_domain):
        version = 1
        for component in (1, 2, 3):
            # Map nodal values + xi_3 derivative to same value,
            if hasDerivatives:
                mappedDerivatives = [
                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
            else:
                mappedDerivatives = [
                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV]
            for derivative in mappedDerivatives:
                boundaryConditions.ConstrainNodeDofsEqual(dependentField, iron.FieldVariableTypes.U,
                        version, derivative, component, nodes, 1.0)
            # Fix xi_1 derivative to be zero
            if hasDerivatives:
                for node in nodes:
                    for derivative in [
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                            iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,
                            ]:
                        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
                                version, derivative, node, component,
                                iron.BoundaryConditionsTypes.FIXED, 0.0)
    elif any(in_domain):
        raise RuntimeError("Mesh decomposition has DOFs "
                "that must be constrained to be equal in separate domains")

# Apply incremented cavity pressure on the endocardial surface:
internalNodes = set(geometry.nodeGroup('internal'))
for node in internalNodes.intersection(geometricDomainNodes):
    derivative, version = 1, 1
    # xi_3 is the transmural direction
    xiDirection = 3
    # For pressure/force boundary conditions, the DELUDELN field variable is
    # constrained rather than the U field variable
    boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.DELUDELN,
            version, derivative, node, xiDirection,
            iron.BoundaryConditionsTypes.PRESSURE_INCREMENTED, -cavityPressure)

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Copy deformed geometry into deformed field
for component in [1, 2, 3]:
    dependentField.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component,
        deformedField, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component)

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results to exnode/exelem files
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/prolate_spheroid", "FORTRAN")
fields.ElementsExport("./results/prolate_spheroid", "FORTRAN")
fields.Finalise()

# Finalise OpenCMISS-Iron
iron.Finalise()
