from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
meshPath = "/home/lmmp/thixotropicSimulation/PreProcessing/ParallelPlates/"
meshFile = "ParallelPlates"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
# mesh.msh2hdmf3D()
print(mesh.subdomains)
mesh.createMeshObject3D()

boundaries = Boundaries(mesh=mesh, Pin=1)

fluid = Fluid.from_file("/home/lmmp/thixotropicSimulation/PreProcessing/ParallelPlates/fluid.txt")

newtonianTest = Solver(mesh=mesh,fluid=fluid,boundaries=boundaries)
newtonianTest.NewtonianSolver()
# newtonianTest.SaveSimulationData("ParallelPlates3D")
