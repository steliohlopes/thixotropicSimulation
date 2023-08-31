from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
meshPath = "/home/lmmp/thixotropicSimulation/PreProcessing/PipeFlow2D/"
meshFile = "PipeFlow2D"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
# mesh.msh2hdmf3D()
print(mesh.subdomains)
mesh.createMeshObject2D()

boundaries = Boundaries(mesh=mesh, Pin=1)

fluid = Fluid(
        rho=1000,
        k=0.1,
        nPow=0.7,
        tau0=6.21358,
        eta0=0.001,
        etaInf=64.1,
        ts=663
        )
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)
# problem.NewtonianEquation()
# newtonianTest = Solver(problem)
# newtonianTest.SimulateEquation()

# problem.PowerLawEquation()
# newtonianTest.SimulateEquation()
problem.ThixotropicEquation()
tixotropicTest = Solver(problem)
# newtonianTest.SaveSimulationData(filePath=meshPath,fileName="ParallelPlatesNewtonian")

