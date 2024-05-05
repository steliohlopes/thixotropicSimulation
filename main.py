from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/ParallelPlatesSymmetry/"
meshFile = "ParallelPlatesSymmetry"
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

fluid = Fluid(
        rho=1000,
        k=1,
        nPow=0.8,
        phi0=0.001,
        phiInf=2,
        Ta = 1,
        Tc = 1
        )

Pin = 2000
Pout = 0
L = 2
H=0.2
U = -((Pout-Pin)/L)*(H/2/(2*fluid.k))*(H-H/2)

boundaries = Boundaries(mesh=mesh, Pin=Pin,Pout=Pout,symmetryBCs=["Symmetry"],symmetryAxis=2)

problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)

problem.GNFEquation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="ParallelPlatesSymmetryNewtonian")
wini = problem.w
del newtonianTest

problem.GNFEquation('SMD')
PowerLawTest = Solver(problem,absTol = 1e-6)
PowerLawTest.SimulateEquation()
PowerLawTest.SaveSimulationData(filePath=meshPath,fileName="ParallelPlatesSymmetrySMD")

# boundaries.change_parameter(Fluidityin=0.1)
# problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)
# problem2.ThixotropicEquation(wini=wini)
# del wini
# Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
# Thixotropic.SimulateEquation()
# Thixotropic.SaveSimulationData(filePath=meshPath,fileName="CoatingBarSymmetrytinyThixotropicLaponite")
# wini = problem2.w
# del Thixotropic, problem2

# fluid2 = Fluid(
#         rho=1000,
#         k=1,
#         nPow=0.32,
#         phi0=0.001,
#         phiInf=50,
#         Ta = 10,
#         Tc = 10
#         )


# problem3 = Problem(mesh=mesh,fluid=fluid2,boundaries=boundaries)

# problem3.ThixotropicEquation(wini=wini)
# del wini
# Thixotropic2 = Solver(problem3,maxIter = 100,absTol = 1e-6)
# Thixotropic2.SimulateEquation()
# Thixotropic2.SaveSimulationData(filePath=meshPath,fileName="CoatingBarSymmetrytinyThixotropicLaponite")