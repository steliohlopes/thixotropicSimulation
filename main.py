from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/PipeFlowPlate3D/"
meshFile = "PipeFlowPlate3D"

mesh = FiniteElementMesh(meshPath=meshPath,meshFile=meshFile)
if comm.rank ==0:
        print(mesh.subdomains)
        if not os.path.exists(f'{meshPath}mesh.xdmf'):
                mesh.msh2hdmf3D()
                
mesh.createMeshObject3D()

info("Num DOFs {}".format(mesh.DoF))

boundaries = Boundaries(mesh=mesh, Pin=1e5)

fluid = Fluid(
        rho=1000,
        k=1,
        nPow=0.60,
        phi0=0.001,
        phiInf=60,
        Ta = 1,
        Tc = 1
        )
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries)

problem.GNFEquation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="CoatingBarNewtonian")
# newtonianTest.viscoty_plot(arrayx=[-4e-3,4e-3],fileName=f'{meshPath}NewtonianViscosity.png')

problem.GNFEquation('SMD')
PowerLawTest = Solver(problem,maxIter = 100)
PowerLawTest.SimulateEquation()
PowerLawTest.SaveSimulationData(filePath=meshPath,fileName="CoatingBarSMD")
# PowerLawTest.velocity_plot(R=100e-6,xpoint= 23e-3/2,fileName=f'{meshPath}SMDResult.png')
# PowerLawTest.viscoty_plot(arrayx=[-4e-3,4e-3],fileName=f'{meshPath}GNFViscosity.png')

boundaries.change_parameter(Fluidityin=0.1)
problem.ThixotropicEquation()
Thixotropic = Solver(problem,maxIter = 100)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName="CoatingBarThixotropic")
# Thixotropic.velocity_plot(R=100e-6,xpoint= 23e-3/2,fileName=f'{meshPath}ThixotropicResultTa{problem.fluid.Ta}.png')
# Thixotropic.viscoty_plot(arrayx=[-4e-3,4e-3],fileName=f'{meshPath}ThixotropicViscosity{problem.fluid.Ta}.png')