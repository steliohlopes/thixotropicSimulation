from PreProcessing.mesh import FiniteElementMesh
from ProblemInputs import Fluid
from Solver.boundaries import Boundaries
from Solver.equations import Solver
from Solver.problem import Problem
from dolfin import *
import os

comm = MPI.comm_world

meshPath = "/home/stelio/thixotropicSimulation/PreProcessing/CoatingHangerMeng/"
meshFile = "CoatingHangerMeng"
simulation_type = '3D'

R = 10e-3
L=168e-3
H=1.5e-3

Q0 = 0.00011371668940800001 #[m3/s]

Q = Q0

#TODO CALCULAR DE NOVO
meanArea=1.8326e-4

Lambda=1

Origin=[-0.07,0,0]

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
        phi0=1/50,
        phiInf=1/0.01,
        Ta = Lambda*meanArea*L/Q,
        Tc =  Lambda*meanArea*L/Q
        )


####################################################################################
# NEWTONIAN
####################################################################################
# checkPointFile = f"{meshPath}{meshFile}Newtonian_checkpoint.xdmf"
boundaries = Boundaries(mesh=mesh,fluid=fluid, Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries, L=H,Q=Q)
problem.Equation('newtonian')
newtonianTest = Solver(problem)
newtonianTest.SimulateEquation()
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Newtonian",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}SMD3",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}SMD4",dimensional=True,Checkpoint=True)

newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic5a",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic5b",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic5c",dimensional=True,Checkpoint=True)

newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic6a",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic6b",dimensional=True,Checkpoint=True)
newtonianTest.SaveSimulationData(filePath=meshPath,fileName=f"{meshFile}Thixotropic6c",dimensional=True,Checkpoint=True)

del problem,newtonianTest


######################################################################################
#SMD
######################################################################################
Q=Q0*5
boundaries = Boundaries(mesh=mesh,fluid=fluid, Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries, L=H,Q=Q)

checkPointFile = f"{meshPath}{meshFile}SMD3_checkpoint.xdmf"
caseName="SMD3"
minNpow_achieved = 1 # Track minimum successfully simulated nPow
objectiveNpow = 0.696

problem.fluid.nPow = 0.696
maxIter = 10
iter= 0
while problem.fluid.nPow >= objectiveNpow and iter<=maxIter :
    try:
        iter+=1
        if comm.rank ==0:
                print("=" * 100) 
                info("case {}".format(caseName))
                info("NPow {}".format(fluid.nPow))
        problem.Equation('SMD', CheckPointFile=checkPointFile)
        newtonianTest = Solver(problem)
        newtonianTest.SimulateEquation()
        newtonianTest.SaveSimulationData(filePath=meshPath, fileName=f"{meshFile}{caseName}", dimensional=True, Checkpoint=True)
        newtonianTest.velocity_plot(sweep_dict=sweep_dict, velocity_coord=velocity_coord, num_points=num_points, fileName=f"{meshPath}{meshFile}{caseName}.csv")
        if problem.fluid.nPow == objectiveNpow:
              break
        # Update minNpow_achieved only if successful
        minNpow_achieved = problem.fluid.nPow
        iter=0
        # Controlled adjustment of nPow (replace with your desired logic)
        adjustment = min(0.11, problem.fluid.nPow - objectiveNpow)  # Adjust by max 0.1
        problem.fluid.nPow -= adjustment

    except:
        if problem.fluid.nPow <= minNpow_achieved:  # Use the achieved value
            problem.fluid.nPow = (minNpow_achieved + problem.fluid.nPow) / 2
        else:
            break

del problem,newtonianTest, boundaries


######################################################################################
#SMD
######################################################################################
Q=Q0/5
boundaries = Boundaries(mesh=mesh,fluid=fluid, Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])
problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries, L=H,Q=Q)

checkPointFile = f"{meshPath}{meshFile}SMD4_checkpoint.xdmf"
caseName="SMD4"
minNpow_achieved = 1  # Track minimum successfully simulated nPow
objectiveNpow = 0.696

problem.fluid.nPow = 0.9
maxIter = 10
iter= 0
while problem.fluid.nPow >= objectiveNpow and iter<=maxIter :
    try:
        iter+=1
        if comm.rank ==0:
                print("=" * 100) 
                info("case {}".format(caseName))
                info("NPow {}".format(fluid.nPow)) 
        problem.Equation('SMD', CheckPointFile=checkPointFile)
        newtonianTest = Solver(problem)
        newtonianTest.SimulateEquation()
        newtonianTest.SaveSimulationData(filePath=meshPath, fileName=f"{meshFile}{caseName}", dimensional=True, Checkpoint=True)
        newtonianTest.velocity_plot(sweep_dict=sweep_dict, velocity_coord=velocity_coord, num_points=num_points, fileName=f"{meshPath}{meshFile}{caseName}.csv")
        if problem.fluid.nPow == objectiveNpow:
              break
        # Update minNpow_achieved only if successful
        minNpow_achieved = problem.fluid.nPow
        iter=0
        # Controlled adjustment of nPow (replace with your desired logic)
        adjustment = min(0.11, problem.fluid.nPow - objectiveNpow)  # Adjust by max 0.1
        problem.fluid.nPow -= adjustment

    except:
        if problem.fluid.nPow <= minNpow_achieved:  # Use the achieved value
            problem.fluid.nPow = (minNpow_achieved + problem.fluid.nPow) / 2
        else:
            break

del problem,newtonianTest

######################################################################################
#THIXOTROPIC
######################################################################################
Q=Q0/5
caseName5 = ["Thixotropic5a","Thixotropic5b","Thixotropic5c"]
LambdaOpt5=[5e-2,1,5]

for i in range(len(caseName5)):
        case=caseName5[i]
        Lambda=LambdaOpt5[i]

        boundaries = Boundaries(mesh=mesh,fluid=fluid, Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])
        problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries, L=H,Q=Q)
        checkPointFile = f"{meshPath}{meshFile}{case}_checkpoint.xdmf"
        problem.fluid.Ta = Lambda*meanArea*L/Q
        problem.fluid.Tc = Lambda*meanArea*L/Q
        minNpow_achieved = 1 # Track minimum successfully simulated nPow
        objectiveNpow = 0.696

        problem.fluid.nPow = 0.9
        maxIter = 10
        iter= 0
        while problem.fluid.nPow >= objectiveNpow and iter<=maxIter :
                try:
                        if comm.rank ==0:
                                print("=" * 100) 
                                info("case {}".format(case))
                                info("NPow {}".format(fluid.nPow)) 
                        iter+=1
                        problem.Equation('thixotropic', CheckPointFile=checkPointFile)
                        newtonianTest = Solver(problem)
                        newtonianTest.SimulateEquation()
                        newtonianTest.SaveSimulationData(filePath=meshPath, fileName=f"{meshFile}{case}", dimensional=True, Checkpoint=True)
                        newtonianTest.velocity_plot(sweep_dict=sweep_dict, velocity_coord=velocity_coord, num_points=num_points, fileName=f"{meshPath}{meshFile}{case}.csv")
                        if problem.fluid.nPow == objectiveNpow:
                                break
                        # Update minNpow_achieved only if successful
                        minNpow_achieved = problem.fluid.nPow
                        iter=0
                        # Controlled adjustment of nPow (replace with your desired logic)
                        adjustment = min(0.11, problem.fluid.nPow - objectiveNpow)  # Adjust by max 0.1
                        problem.fluid.nPow -= adjustment

                except:
                        if problem.fluid.nPow <= minNpow_achieved:  # Use the achieved value
                                problem.fluid.nPow = (minNpow_achieved + problem.fluid.nPow) / 2
                        else:
                                break

        del problem,newtonianTest,boundaries
######################################################################################
#THIXOTROPIC 6
######################################################################################

Q=Q0*5
caseName6 = ["Thixotropic6a","Thixotropic6b","Thixotropic6c"]
LambdaOpt6=[5e-2,1,5]

for i in range(len(caseName6)):
        case=caseName6[i]
        Lambda=LambdaOpt6[i]

        boundaries = Boundaries(mesh=mesh,fluid=fluid, Origin=Origin,R=R,symmetryBCs=[['SymmetryY'],['SymmetryZ']],symmetryAxis=[1,2])
        problem = Problem(mesh=mesh,fluid=fluid,boundaries=boundaries, L=H,Q=Q)
        checkPointFile = f"{meshPath}{meshFile}{case}_checkpoint.xdmf"
        problem.fluid.Ta = Lambda*meanArea*L/Q
        problem.fluid.Tc = Lambda*meanArea*L/Q
        minNpow_achieved = 1 # Track minimum successfully simulated nPow
        objectiveNpow = 0.696

        problem.fluid.nPow = 0.9
        maxIter = 10
        iter= 0
        while problem.fluid.nPow >= objectiveNpow and iter<=maxIter :
                try:
                        if comm.rank ==0:
                                print("=" * 100) 
                                info("case {}".format(case))
                                info("NPow {}".format(fluid.nPow)) 
                        iter+=1
                        problem.Equation('thixotropic', CheckPointFile=checkPointFile)
                        newtonianTest = Solver(problem)
                        newtonianTest.SimulateEquation()
                        newtonianTest.SaveSimulationData(filePath=meshPath, fileName=f"{meshFile}{case}", dimensional=True, Checkpoint=True)
                        newtonianTest.velocity_plot(sweep_dict=sweep_dict, velocity_coord=velocity_coord, num_points=num_points, fileName=f"{meshPath}{meshFile}{case}.csv")
                        if problem.fluid.nPow == objectiveNpow:
                                break
                        # Update minNpow_achieved only if successful
                        minNpow_achieved = problem.fluid.nPow
                        iter=0
                        # Controlled adjustment of nPow (replace with your desired logic)
                        adjustment = min(0.11, problem.fluid.nPow - objectiveNpow)  # Adjust by max 0.1
                        problem.fluid.nPow -= adjustment

                except:
                        if problem.fluid.nPow <= minNpow_achieved:  # Use the achieved value
                                problem.fluid.nPow = (minNpow_achieved + problem.fluid.nPow) / 2
                        else:
                                break

        del problem,newtonianTest,boundaries