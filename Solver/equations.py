import sys

sys.path.append("..")
from dolfin import *
from ProblemInputs import Inputs
from PreProcessing.mesh import FiniteElementMesh
from Solver.boundaries import Boundaries

# Deformation Tensor
def DD(u):
    # Cartesian
    D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    return D


# Stress Tensor
def TT(u, p, mu):
    # Cartesian
    T = 2 * mu * DD(u) - p * Identity(len(u))
    return T


def gammaDot(u):
    return pow(2 * inner(DD(u), DD(u)), 0.5)


def eta(k,nPow,u):
    eps=DOLFIN_EPS
    return k*pow(gammaDot(u)+eps,nPow-1)

class Solver:
    def __init__(self, input=Inputs(), mesh=FiniteElementMesh(), boundaries = Boundaries()):
        
        #TODO HERDAR COM SUPER
        self.k = input.k
        self.Pout = input.Pout
        self.outletBCs = input.outletBCs

        ## Solver Parameters
        self.absTol = input.absTol
        self.relTol = input.relTol
        self.maxIter = input.maxIter
        self.linearSolver = input.linearSolver
        self.nonlinearSolver = input.nonlinearSolver

        self.subdomains = mesh.subdomains
        self.dx = mesh.dx
        self.ds = mesh.ds
        self.n = mesh.n

        self.bcs = boundaries.bcs
        self.inletCondition = boundaries.inletCondition
        self.inletBCs = boundaries.inletBCs

        if self.inletCondition == 0:
            self.Pin == boundaries.Pin
        elif self.inletCondition == 1:
            self.Uin == boundaries.Uin

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction(mesh.functionSpace)
        (self.v, self.q) = TestFunctions(mesh.functionSpace) #Funçao peso
        self.w = Function(mesh.functionSpace)

        # Split into Velocity and Pressure
        if mesh.Dim == 3:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1], self.w[2])), self.w[3])
        elif mesh.Dim == 2:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1])), self.w[2])


    def NewtonianSolver(self):
        a01 = (inner(TT(self.u,self.p,eta(self.k,1,self.u)),DD(self.v)))*self.dx()
        # + (rho*dot(dot(u,grad(u)),v) 
                                     
        L01 =  - (self.Pout)*dot(self.n,self.v)*self.ds(self.subdomains[self.outletBCs[0]]) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.inletCondition == 0:
            L01 = L01 - (self.Pin)*dot(self.n,self.v)*self.ds(self.subdomains[self.inletBCs[0]]) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.dx()
        L02 = 0

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        problemU0 = NonlinearVariationalProblem(F0,self.w,self.bcs,J0)
        solverU0 = NonlinearVariationalSolver(problemU0)
            # # Solver Parameters
        prmU0 = solverU0.parameters 
        prmU0['nonlinear_solver'] = self.nonlinearSolver
        prmU0['newton_solver']['absolute_tolerance'] = self.absTol
        prmU0['newton_solver']['relative_tolerance'] = self.relTol
        prmU0['newton_solver']['maximum_iterations'] = self.maxIter
        prmU0['newton_solver']['linear_solver'] = self.linearSolver
        prmU0['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        (no_iterations,converged) = solverU0.solve()

        return self.w
    




        
    


