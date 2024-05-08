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

Pin = 3000
Pout = 0
W_outlet=10e-3
L=W_outlet/2
# H=100e-6*2
# U = -((Pout-Pin)/L)*(H/2/(2*fluid.k))*(H-H/2)
U=1e-3

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

boundaries = Boundaries(mesh=mesh, Pin=Pin,Pout=Pout,symmetryBCs=["Symmetry"],symmetryAxis=2)

problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)

problem.Equation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName="CoatingHangerSymmetry2Newtonian",dimensional=True)
newtonianTest.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName="CoatingHangerSymmetry2Newtonian.png")
wini = problem.w
del newtonianTest

boundaries.change_parameter(Fluidityin=0.1)
problem2 = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
problem2.Equation(wini=wini,model='thixotropic')

times = [5,1,0.1,0.01,0.001,0.0001]

Thixotropic = Solver(problem2,maxIter = 100,absTol = 1e-6)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName=f"CoatingHangerSymmetry2Thixotropic{problem2.fluid.Ta}",dimensional=True)
Thixotropic.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetry2Thixotropic{problem2.fluid.Ta}.png")
wini = problem2.w
del Thixotropic

for t in times:
        fluid.Ta=t
        fluid.Tc=t
        problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
        problem.Equation(wini=wini,model='thixotropic')
        Thixotropic = Solver(problem,maxIter = 100,absTol = 1e-6)
        try:
                Thixotropic.SimulateEquation()
                Thixotropic.SaveSimulationData(filePath=meshPath,fileName=f"CoatingHangerSymmetry2Thixotropic{t}",dimensional=True)
                Thixotropic.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName=f"CoatingHangerSymmetry2Thixotropic{t}.png")
                wini = problem.w
                del Thixotropic
        except:
                pass

del boundaries
boundaries = Boundaries(mesh=mesh, Pin=Pin,Pout=Pout,symmetryBCs=["Symmetry"],symmetryAxis=2)  
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries,U = U, L=L)
problem.Equation(wini=wini,model='SMD')
Thixotropic = Solver(problem,maxIter = 100,absTol = 1e-6)
Thixotropic.SimulateEquation()
Thixotropic.SaveSimulationData(filePath=meshPath,fileName=f"CoatingHangerSymmetry2SMD",dimensional=True)
Thixotropic.velocity_plot(sweep_dict=sweep_dict,velocity_coord=velocity_coord,num_points=num_points,fileName="CoatingHangerSymmetry2SMD.png")