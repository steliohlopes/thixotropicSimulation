from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/CoatingHangerSymmetry2/"
meshFile = "CoatingHangerSymmetry2"
simulation_type = '3D'

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)

if comm.rank ==0:
        print(mesh.subdomains)
        if not os.path.exists(f'{meshPath}mesh.xdmf'):
                if simulation_type =="3D":
                        mesh.msh2hdmf3D()
                else:
                        mesh.msh2hdmf2D()

if simulation_type =="3D":
        mesh.createMeshObject3D()
else:
        mesh.createMeshObject2D()
        
if comm.rank ==0:
        info("Num DOFs {}".format(mesh.DoF))         

boundaries = Boundaries(mesh=mesh, Pin=2000,symmetryBCs=["Symmetry"],symmetryAxis=2)

fluid = Fluid(
        rho=1000,
        k=1,
        nPow=0.9,
        phi0=0.001,
        phiInf=50,
        Ta = 10,
        Tc = 10
        )
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)

problem.GNFEquation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="CoatingHangerSymmetry2Newtonian")
wini = problem.w
del newtonianTest,problem

# problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)
# problem.GNFEquation('SMD',wini=wini)
# PowerLawTest = Solver(problem,maxIter = 100,absTol = 1e-9)
# PowerLawTest.SimulateEquation()
# PowerLawTest.SaveSimulationData(filePath=meshPath,fileName="CoatingHangerSymmetry2SMD")

boundaries.change_parameter(Fluidityin=0.1)
problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)
problem2.ThixotropicEquation(wini=wini)
del wini
Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName="CoatingHangerSymmetry2Thixotropic3")
wini = problem2.w
del Thixotropic, problem2

fluid2 = Fluid(
        rho=1000,
        k=1,
        nPow=0.85,
        phi0=0.001,
        phiInf=50,
        Ta = 2,
        Tc = 2
        )


problem3 = Problem(mesh=mesh,fluid=fluid2,boundaries=boundaries)

problem3.ThixotropicEquation(wini=wini)
del wini
Thixotropic2 = Solver(problem3,maxIter = 100,absTol = 1e-6)
Thixotropic2.SimulateEquation()
Thixotropic2.SaveSimulationData(filePath=meshPath,fileName="CoatingHangerSymmetry2Thixotropic2")