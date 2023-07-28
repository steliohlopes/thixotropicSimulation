import sys

sys.path.append("..")
from dolfin import *
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
    def __init__(self, mesh,fluid, boundaries):

        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries

        ## Solver Parameters
        self.nonlinearSolver = "newton"
        self.absTol = 1e-9
        self.relTol = 1e-10
        self.maxIter = 30
        self.linearSolver = "mumps"

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction( self.mesh.functionSpace)
        (self.v, self.q) = TestFunctions( self.mesh.functionSpace) #Funçao peso
        self.w = Function( self.mesh.functionSpace)

        # Split into Velocity and Pressure
        if  self.mesh.Dim == 3:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1], self.w[2])), self.w[3])
        elif  self.mesh.Dim == 2:
            (self.u, self.p) = (as_vector((self.w[0], self.w[1])), self.w[2])


    def NewtonianSolver(self,wini = None):
        if wini != None:
            self.w = wini

        a01 = (inner(TT(self.u,self.p,eta(self.fluid.k,1,self.u)),DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 
                                     
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.outletBCs[0]]) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.inletBCs[0]]) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)
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
    
    def PowerLawSolver(self,wini = None):
        if wini != None:
            self.w = wini
        
        a01 = (inner(TT(self.u,self.p,eta(self.fluid.k,self.fluid.nPow,self.u)),DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 
                                     
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.outletBCs[0]]) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(self.mesh.subdomains[self.boundaries.inletBCs[0]]) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)
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




        
    


