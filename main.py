from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/CoatingHangerSymmetry/"
meshFile = "CoatingHangerSymmetry"
simulation_type = '3D'

# L = 24e-3
# U=1e-5
# Origin=[-L/2,0]
# R = 100e-6


W_outlet=10e-3
L=W_outlet/2
U=5e-2
D_inlet = 0.5e-3
R = D_inlet/2
L_inlet = 1e-3

Origin=[-L_inlet,-R,0]


sweep_dict={0:0.005,1:-0.0001,2:[0,0.005]}
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
        k=1,
        nPow=0.8,
        phi0=0.001,
        phiInf=20,
        Ta = 10,
        Tc = 10
        )


if comm.rank ==0:
        info("L characteristic {}".format(L))
        info("U characteristic {}".format(U))

boundaries = Boundaries(mesh=mesh, UinMax_dim=1,Origin=Origin,R=R,symmetryBCs=['Symmetry'],symmetryAxis=2)

problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)

problem.Equation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"CoatingHangerSymmetryNewtonian",dimensional=True)
newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetryThixotropic{problem.fluid.Ta}.csv")

wini = problem.w
del problem, newtonianTest

times = [10,5,1,0.1]
for i in times:
        fluid.Ta=i
        fluid.Tc=i
        try:
                problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
                problem.Equation('thixotropic',wini=wini)
                newtonianTest = Solver(problem)
                newtonianTest.SimulateEquation()
                newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"CoatingHangerSymmetryThixotropic{problem.fluid.Ta}",dimensional=True)
                newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetryThixotropic{problem.fluid.Ta}.csv")
                wini=problem.wini
        except:
                continue



# boundaries.change_parameter(Fluidityin=0.1)
# problem.Equation('SMD')
# newtonianTest = Solver(problem)
# newtonianTest.SimulateEquation()
# newtonianTest.SaveSimulationData(filePath=meshPath,fileName="ConstrictedSMD",dimensional=False)



# problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L,Pinf=Pinf)
# problem2.Equation(wini=wini,model='thixotropic')

# Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
# Thixotropic.SimulateEquation()
# Thixotropic.SaveSimulationData(filePath=meshPath,fileName=f"ConstrictedThixotropic{problem2.fluid.Ta}",dimensional=True)
# # Thixotropic.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetry2Thixotropic{problem2.fluid.Ta}.png")
