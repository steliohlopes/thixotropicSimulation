from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/CoatingBarSymmetry/"
meshFile = "CoatingBarSymmetry"
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

boundaries = Boundaries(mesh=mesh, Pin=1e5,symmetryBCs=["Symmetry"],symmetryAxis=2)

fluid = Fluid(
        rho=1000,
        k=1,
        nPow=0.60,
        phi0=0.001,
        phiInf=15,
        Ta = 1,
        Tc = 1
        )
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)

problem.GNFEquation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="CoatingBarSymmetryNewtonian")


# problem.GNFEquation('SMD')
# PowerLawTest = Solver(problem,maxIter = 100)
# PowerLawTest.SimulateEquation()
# PowerLawTest.SaveSimulationData(filePath=meshPath,fileName="PipeFlowSMD")


boundaries.change_parameter(Fluidityin=0.1)
problem.ThixotropicEquation()
Thixotropic = Solver(problem,maxIter = 100)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName="CoatingBarSymmetryThixotropic")
