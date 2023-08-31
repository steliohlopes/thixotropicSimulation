from dolfin import *
import sys
sys.path.append("..")
from math import exp



class Problem:
    def __init__(self, mesh,fluid, boundaries):
        self.mesh = mesh
        self.fluid = fluid
        self.boundaries = boundaries

        ##### Functions
        ## Trial and Test function(s)
        self.dw = TrialFunction(self.mesh.functionSpace)
        (self.v, self.q,self.m) = TestFunctions(self.mesh.functionSpace)
        self.w = Function(self.mesh.functionSpace)

        # Split into Velocity and Pressure
        if  self.mesh.Dim == 3:
            (self.u, self.p, self.f) = (as_vector((self.w[0], self.w[1], self.w[2])), self.w[3], self.w[4])
        elif  self.mesh.Dim == 2:
            (self.u, self.p, self.f) = (as_vector((self.w[0], self.w[1])), self.w[2], self.w[3])
    
    # Deformation Tensor
    def DD(self, u):
        # Cartesian
        D = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
        return D


    # Stress Tensor
    def TT(self, u, p, mu):
        # Cartesian
        T = 2 * mu * self.DD(u) - p * Identity(len(u))
        return T


    def gammaDot(self, u):
        return pow(2 * inner(self.DD(u), self.DD(u)), 0.5)


    def eta(self, k,nPow,u):
        return k*pow(self.gammaDot(u)+DOLFIN_EPS,nPow-1)
    
    def sigmoid(self,field):
        # a=50000
        f_expr = Expression("1/(1+exp(-5000*field))", field=field, degree=1)
        # H = 1/(1+exp(-a*x))
        f = interpolate(f_expr, field.function_space())
        return f
    
    def phieq(self,k,nPow,phi0,phiInf,u,p,phiLocal):
        sigma = pow(0.5*tr((self.TT(u,p,(1/phiLocal))*(self.TT(u,p,(1/phiLocal))).T)),0.5)
        b = (pow(sigma/k,(1/nPow))/sigma)
        PHIeq = b/((phiInf-phi0)+b)
        return PHIeq
    
    def Tc(self):
        tc = 663
        return tc
    
    def Ta(self):
        ta = 200
        return ta
    
    def S(self,phieq):
        s = (8/(exp(phieq/0.09)-1))+1.2
        return s


    def NewtonianEquation(self,wini = None):
        if wini != None:
            self.w = wini

        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,1,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)    
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Mass Conservation(Continuity)
        a03 = self.f
        L03 = 1/self.fluid.k

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
    
    def PowerLawEquation(self,wini = None):
        if wini != None:
            self.w = wini
        
        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,self.fluid.nPow,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)                             
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        # Complete Weak Form
        F0 = (a01 + a02) - (L01 + L02)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
    
    def ThixotropicEquation(self,wini = None):
        if wini != None:
            self.w = wini

        a01 = (inner(self.TT(self.u,self.p,self.eta(self.fluid.k,self.fluid.nPow,self.u)),self.DD(self.v)))*self.mesh.dx()
        # + (rho*dot(dot(u,grad(u)),v) 

        outletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.outletBCs if key in self.mesh.subdomains)                             
        L01 =  - (self.boundaries.Pout)*dot(self.mesh.n,self.v)*self.mesh.ds(outletBCsIndex) # Outlet Pressure
            # + inner(rho*fb(inputs),v)*dx()   # Gravity

        if self.boundaries.inletCondition == 0:
            inletBCsIndex = tuple(self.mesh.subdomains[key] for key in self.boundaries.inletBCs if key in self.mesh.subdomains)
            L01 = L01 - (self.boundaries.Pin)*dot(self.mesh.n,self.v)*self.mesh.ds(inletBCsIndex) # Inlet Pressure 

        # Mass Conservation(Continuity)
        a02 = (self.q*div(self.u))*self.mesh.dx()
        L02 = 0

        #Fluidity
        a03 = (
            inner(self.u,grad(self.f))+
            inner((self.f-self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))/self.Tc(),
             self.sigmoid(self.f-self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)))  -
            (1-self.sigmoid(self.f-self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)))*
            self.S(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))/(self.Ta()/self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))*
            pow((self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)-self.f),(self.S(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))+1)/self.S(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)))*
            pow(self.f,(self.S(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f))-1)/self.S(self.phieq(self.fluid.k,self.fluid.nPow,self.fluid.phi0,self.fluid.phiInf,self.u,self.p,self.f)))
               )*self.m*self.mesh.dx()
        L03 = 0 

        # Complete Weak Form
        F0 = (a01 + a02 + a03) - (L01 + L02 + L03)
        # Jacobian Matrix
        J0 = derivative(F0,self.w,self.dw)

        # Problem and Solver definitions
        self.problemU0 = NonlinearVariationalProblem(F0,self.w,self.boundaries.bcs,J0)

        return self.problemU0
