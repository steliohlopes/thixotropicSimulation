from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
meshPath = "/home/lmmp/thixotropicSimulation/PreProcessing/Contraction/"
meshFile = "Contraction"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
# mesh.msh2hdmf2D()
print(mesh.subdomains)
mesh.createMeshObject2D()

boundaries = Boundaries(mesh=mesh, Pin=1e6,Fluidityin=0.5)

fluid = Fluid(
        rho=1000,
        k=1,
        nPow=0.65,
        tau0=6.21358,
        eta0=0.001,
        etaInf=64.1,
        ts=663,
        phi0=0.001,
        phiInf=15
        )
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)

problem.NewtonianEquation()
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="ContractionNewtonian")

# problem.PowerLawEquation()
# PowerLawTest = Solver(problem,maxIter = 100)
# PowerLawTest.SimulateEquation()
# PowerLawTest.SaveSimulationData(filePath=meshPath,fileName="ContractionPowerLaw")


problem.ThixotropicEquation()
tixotropicTest = Solver(problem,absTol=1.5e-12,maxIter = 100)
tixotropicTest.SimulateEquation()
tixotropicTest.SaveSimulationData(filePath=meshPath,fileName="ContractionThixotropic")

# problem.fluid.nPow=0.30
# problem.ThixotropicEquation()
# tixotropicTest = Solver(problem,absTol=1.5e-12,maxIter = 100)
# tixotropicTest.SimulateEquation()
# tixotropicTest.SaveSimulationData(filePath=meshPath,fileName="ContractionThixotropic")


