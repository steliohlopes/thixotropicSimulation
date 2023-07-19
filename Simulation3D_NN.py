# Libraries Import
import matplotlib.pyplot as plt
import numpy as np
from dolfin import *
import os
import timeit
comm = MPI.comm_world
rank = MPI.rank(comm)

# Helping self defined Functions
# Read subdomains from .msh file
def readDomains(inPath,inFile):
    # Read .msh File
    fid = open(inPath+inFile+'.msh', 'r')
    # Initialize variables
    found = 0
    finished = 0
    physicalNames = {}
    # Loop througn .msh lines
    for line in fid:
        if '$EndPhysicalNames' in line:
            finished == 1
            break
        elif '$PhysicalNames' in line:
            found = 1
        elif found==1 and finished == 0:
            word=line.split()
            if len(word)==3:
                physicalNames[word[2][1:len(word[2])-1]] = int(word[1])

    return physicalNames

# Deformation Tensor
def DD(u):
    #Cartesian
    D = 0.5*(nabla_grad(u) + nabla_grad(u).T)
    return D

# Stress Tensor
def TT(u, p, mu):
    #Cartesian
    T = 2*mu*DD(u) - p*Identity(len(u))
    return T

def gammaDot(u):
    return pow(2*inner(DD(u),DD(u)),0.5)

def eta(k,nPow,u):
    eps=DOLFIN_EPS
    return k*pow(gammaDot(u)+eps,nPow-1)

# Inputs
# Mesh File
meshPath = './Whistle/'
meshFile = 'whistle'

# Pressure Difference or constant velocity
Pin = 0.5
Pout = 0
Uin=0.001

#Geometry


D_inlet=0.5e-3
L_inlet=5e-4
H_inlet=200e-6
H_outlet=H_inlet
W_outlet=2e-2
L_outlet=2e-2
R=1e-3
alpha=45*(pi/180)
l1 = 2*R/(tan(pi/2 - alpha))
l2 = R*cos(alpha)
l3 = R*sin(alpha)



x_outlet=L_inlet+R+l2+l1+L_outlet

inflow_profile = Expression(('Uin*1.5*(1-(( pow(x[1],2) )/( pow(R,2) )))', '0','0'),\
                            degree=2,Uin=Uin,R=(D_inlet/2))




# Fluid Properties
rho = 1000
mu = 0.1

# Mesh Elements
# Velocity
velocityElementfamily = 'Lagrange'
velocityElementOrder = 2
# Pressure
pressureElementfamily = 'Lagrange'
pressureElementOrder = 1

# # Rheology - Modified SMD (Souza Mendes e Dutra (2004)) + Cure(tauY(t)) 
tau0 = 6.21358           # Dinamic Yield Stress               
eta0 = 0.001            # Viscosity Value for Low shear rates
etaInf = 64.1  # Equilibrium Viscosity(Newtonian Plato: Lowgh shear rates)
k = 0.1              # Consistency Index
nPow = 0.7             # Power-law Index
ts = 663             # Caracteristic viscosity buildup time

# Solver Parameters
absTol = 1e-9
relTol = 1e-10
maxIter =30
linearSolver = 'mumps'

Subdomains = readDomains(meshPath,meshFile)

print(Subdomains)
begin("Total running time: %f\n" % (x_outlet))

meshObj = Mesh()
with XDMFFile(meshPath+"mesh.xdmf") as infile:
    infile.read(meshObj)
mvc = MeshValueCollection("size_t", meshObj, 2)
with XDMFFile(meshPath+"mf.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(meshObj, mvc)

mvc2 = MeshValueCollection("size_t", meshObj, 3)
with XDMFFile(meshPath+"mesh.xdmf") as infile:
    infile.read(mvc2, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(meshObj, mvc2)

# Get Element Shape: Triangle, etc...
elementShape = meshObj.ufl_cell()

# Set Mesh Elements
Uel = VectorElement(velocityElementfamily, elementShape, velocityElementOrder) # Velocity vector field
Pel = FiniteElement(pressureElementfamily, elementShape, pressureElementOrder) # Pressure field
Nel = FiniteElement(pressureElementfamily, elementShape, pressureElementOrder) # Viscosity field
UPel = MixedElement([Uel,Pel])

# Define any measure associated with domain and subdomains
dx = Measure('dx', domain=meshObj, subdomain_data=cf)
ds = Measure('ds', domain=meshObj, subdomain_data=mf)

# Vectors Normal to the Mesh
n = FacetNormal(meshObj) # Normal vector to mesh

# Function Spaces: Flow
# Mixed Function Space: Pressure and Velocity
W = FunctionSpace(meshObj,UPel)

##### Functions
## Trial and Test function(s)
dw = TrialFunction(W)
(v, q) = TestFunctions(W)
w = Function(W)

# Split into Velocity and Pressure
(u, p) = (as_vector((w[0], w[1], w[2])), w[3])

# Apply Flow Boundary Conditions
bcU1 = DirichletBC(W.sub(0),Constant((0.0,0.0,0.0)),mf,Subdomains['Wall'])
# bcU3 = DirichletBC(W.sub(0), inflow_profile, mf,Subdomains['Inlet'])
bcs = [bcU1]


#Start timer 
start = timeit.default_timer()
# eta = Constant(mu)

#----------Newtonian Solver------
# a01 = (rho*dot(dot(u,grad(u)),v) + inner(TT(u,p,eta0),DD(v)))*dx()
a01 = (inner(TT(u,p,eta(k,1,u)),DD(v)))*dx()

       # Outlet Pressure                           # Inlet Pressure                           # Gravity
L01 =  - (Pout)*dot(n,v)*ds(Subdomains['Outlet']) - (Pin)*dot(n,v)*ds(Subdomains['Inlet'])# + inner(rho*fb(inputs),v)*dx()    
        
# Mass Conservation(Continuity)
a02 = (q*div(u))*dx()
L02 = 0

# Complete Weak Form
F0 = (a01 + a02) - (L01 + L02)
    # Jacobian Matrix
J0 = derivative(F0,w,dw)

        
##########   Numerical Solver Properties
# Problem and Solver definitions
problemU0 = NonlinearVariationalProblem(F0,w,bcs,J0)
solverU0 = NonlinearVariationalSolver(problemU0)
    # # Solver Parameters
prmU0 = solverU0.parameters 
prmU0['nonlinear_solver'] = 'newton'
prmU0['newton_solver']['absolute_tolerance'] = absTol
prmU0['newton_solver']['relative_tolerance'] = relTol
prmU0['newton_solver']['maximum_iterations'] = maxIter
prmU0['newton_solver']['linear_solver'] = linearSolver
prmU0['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
# prmU0['newton_solver']['preconditioner'] = 'ilu'
# info(prmU0,True)  #get full info on the parameters
# Solve Problem
(no_iterations,converged) = solverU0.solve()

# # ----------------------------Non Newtonian----------------------------
# wini = interpolate(w, W)
# a01 = (inner(TT(u,p,eta(k,nPow,u)),DD(v)))*dx()
#  # Complete Weak Form
# F0 = (a01 + a02) - (L01 + L02)
#     # Jacobian Matrix
# J0 = derivative(F0,w,dw)
        
# # ##########   Numerical Solver Properties
# # # Problem and Solver definitions
# problemU0 = NonlinearVariationalProblem(F0,w,bcs,J0)
# solverU0 = NonlinearVariationalSolver(problemU0)
# # # # Solver Parameters
# prmU0 = solverU0.parameters 
# prmU0['nonlinear_solver'] = 'newton'
# prmU0['newton_solver']['absolute_tolerance'] = absTol
# prmU0['newton_solver']['relative_tolerance'] = relTol
# prmU0['newton_solver']['maximum_iterations'] = maxIter
# prmU0['newton_solver']['linear_solver'] = linearSolver
# # prmU0['newton_solver']['preconditioner'] = 'ilu'
# # info(prmU0,True)  #get full info on the parameters
# # Solve Problem
# (no_iterations,converged) = solverU0.solve()

(u1, p1) = w.leaf_node().split()
u1.rename("Velocity Vector", "")
p1.rename("Pressure", "")

# End Time
stop = timeit.default_timer()
total_time = stop - start
    
# Output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

begin("Total running time: %dh:%dmin:%ds \n" % (hours, mins, secs))

Simulation_file = XDMFFile(meshPath+meshFile+"NN.xdmf")
Simulation_file.parameters["flush_output"] = True
Simulation_file.parameters["functions_share_mesh"]= True
Simulation_file.write(u1, 0.0)
Simulation_file.write(p1, 0.0)

ux = []
j = []



for i in np.linspace(-0.01, 0.01, 200):
    j.append(i)
    ux.append(u1(0.0242,0.00035,i)[0])

plt.legend(['FEniCS result Velocity Profile'])
plt.show()
plt.savefig(meshPath+'Result_u(y).png')


# #poiseuille
# jPoiseuille=[]
# uxPoiseuille = []

# for i in np.linspace(-R, R, 200):
#     jPoiseuille.append(i)
#     # Poiseuille Pressure bc
#     # uxPoiseuille.append(-((Pout-Pin)/L)*(i/(2*k))*(R-i))
#     #Poiseuille velocity bc y[0,H]
#     # uxPoiseuille.append((1.5*Uin*(1-((i-R/2)**2/(R/2)**2))))
    #Poiseuille velocity bc y[-H/2,H/2]
    # uxPoiseuille.append((1.5*Uin*(1-((i)**2/(R)**2))))

# #PowerLaw Analytical
# jPowerLaw=[]
# uxPowerLaw=[]

# for i in np.linspace(-R, R, 200):
#     jPowerLaw.append(i)
#     uxPowerLaw.append(nPow/(nPow+1)*((p1(xpoint-epsdx,i,0)-p1(xpoint+epsdx,i,0))/(2*k*(2*epsdx)))**(1/nPow)*(R**((nPow+1)/nPow) -abs(i)**((nPow+1)/nPow)))


# plt.plot(
# #     uxPowerLaw, jPowerLaw,'g',
#     uxPoiseuille, jPoiseuille,'b',
#     ux,j,'r')

# plt.legend([
# #     'PowerLaw Analytical',
#     'Poiseuille flow',
#     'FEniCS result'])
# plt.show()
# plt.savefig(meshPath+'Result.png')
