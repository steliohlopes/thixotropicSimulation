from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/Contraction/"
meshFile = "Contraction"
simulation_type = '2D'

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
        phiInf=20,
        Ta = 50,
        Tc = 50
        )

Pin = 2000
Pout = 0
L = 24e-3
H=100e-6*2
U = -((Pout-Pin)/L)*(H/2/(2*fluid.k))*(H-H/2)

if comm.rank ==0:
        info("L characteristic {}".format(L))
        info("U characteristic {}".format(U))

boundaries = Boundaries(mesh=mesh, Pin=Pin,Pout=Pout)

problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)

problem.GNFEquation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="ContractionNewtonian")
wini = problem.w
del newtonianTest


boundaries.change_parameter(Fluidityin=0.01)
problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
problem2.GNFEquation(wini=wini,model='thixotropic')
# del wini
Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName="ContractionThixotropic",dimensional=True)
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