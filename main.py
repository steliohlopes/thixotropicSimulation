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

R = 10e-3
L=168e-3

Q = 0.00011371668940800001 #[m3/s]

U=Q *2/(pi*pow(R,2))

Origin=[-0.01,0,0]

alfa = 30*(pi/180)
slotCenter = L*tan(alfa)

L_outlet = 50e-3
R1 = 10e-3


sweep_dict={0:2*R1+slotCenter+L_outlet,1:0,2:[0,L]}

velocity_coord = 0
num_points = 3000


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
        k=0.799,
        nPow=0.696,
        phi0=0.001,
        phiInf=20,
        Ta = 10,
        Tc = 10
        )


if comm.rank ==0:
        info("L characteristic {}".format(L))
        info("U characteristic {}".format(U))

boundaries = Boundaries(mesh=mesh, UinMax_dim=1,Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])

problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)

# problem.Equation('newtonian')
# newtonianTest = Solver(problem)
# newtonianTest.SimulateEquation()
# newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Newtonian",dimensional=True,Checkpoint=True)
# newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"{meshPath}{meshFile}Newtonian.csv")


# problem.Equation('SMD')
# newtonianTest = Solver(problem)
# newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}SMD",dimensional=True,Checkpoint=True)
# newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"{meshPath}{meshFile}SMD.csv")

# wini = problem.w
# del problem, newtonianTest

times = [10,5,1,0.1]
for i in times:
        fluid.Ta=i
        fluid.Tc=i
        try:
                problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
                problem.Equation('thixotropic',CheckPointFile=f'{meshPath}{meshFile}SMD_checkpoint.xdmf')
                newtonianTest = Solver(problem)
                newtonianTest.SimulateEquation()
                newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic{problem.fluid.Ta}",dimensional=True,Checkpoint=True)
                newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"{meshPath}{meshFile}Thixotropic{problem.fluid.Ta}.csv")
        except:
                continue






# problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L,Pinf=Pinf)
# problem2.Equation(wini=wini,model='thixotropic')

# Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
# Thixotropic.SimulateEquation()
# Thixotropic.SaveSimulationData(filePath=meshPath,fileName=f"ConstrictedThixotropic{problem2.fluid.Ta}",dimensional=True)
# # Thixotropic.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetry2Thixotropic{problem2.fluid.Ta}.png")
